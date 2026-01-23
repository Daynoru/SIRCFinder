// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#ifndef COMPLEXITY_FILTER_H
#define COMPLEXITY_FILTER_H

#include "dna_utils.h"
#include <vector>
#include <string>
#include <mutex>
#include <tbb/task_scheduler_init.h>

using namespace dna_utils;

// Worker 
class ComplexityWorker : public RcppParallel::Worker {
private:
  const std::vector<CoreResult>& cassettes_;
  const std::string& genome_;
  int max_k_;
  double dr_threshold_;
  double spacer_threshold_;
  std::vector<bool>& results_;
  std::mutex& mutex_;
  
public:
  ComplexityWorker(const std::vector<CoreResult>& cassettes,
                   const std::string& genome,
                   int max_k,
                   double dr_threshold,
                   double spacer_threshold,
                   std::vector<bool>& results,
                   std::mutex& mutex)
    : cassettes_(cassettes), genome_(genome), max_k_(max_k),
      dr_threshold_(dr_threshold), spacer_threshold_(spacer_threshold),
      results_(results), mutex_(mutex) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      const auto& cassette = cassettes_[i];
      
      // DR sequences check
      bool dr_passed = true;
      for (const auto& seq : cassette.sequences) {
        if (!passes_complexity_threshold(seq, max_k_, dr_threshold_)) {
          dr_passed = false;
          break;
        }
      }
      
      if (!dr_passed) {
        std::lock_guard<std::mutex> lock(mutex_);
        results_[i] = false;
        continue;
      }
      
      // Spacers check
      auto spacers = extract_spacers(genome_, cassette.positions, cassette.sequences[0].length());
      bool spacer_passed = true;
      
      for (const auto& spacer : spacers) {
        if (!passes_complexity_threshold(spacer, max_k_, spacer_threshold_)) {
          spacer_passed = false;
          break;
        }
      }
      
      std::lock_guard<std::mutex> lock(mutex_);
      results_[i] = dr_passed && spacer_passed;
    }
  }
};



// Main function of complexity estimation
inline std::vector<CoreResult> complexity_filter_internal(
    const std::vector<CoreResult>& cassettes,
    const std::string& genome,
    int max_k,
    double dr_threshold,
    double spacer_threshold,
    int num_threads) {
  
  if (cassettes.empty()) {
    return std::vector<CoreResult>();
  }
  
  // parallel init
  tbb::task_scheduler_init init(num_threads);
  std::vector<bool> results(cassettes.size(), false);
  std::mutex results_mutex;
  
  // worker start
  ComplexityWorker worker(cassettes, genome, max_k, dr_threshold, 
                          spacer_threshold, results, results_mutex);
  
  RcppParallel::parallelFor(0, cassettes.size(), worker);
  
  // Results filtering
  std::vector<CoreResult> filtered_cassettes;
  for (size_t i = 0; i < cassettes.size(); ++i) {
    if (results[i]) {
      filtered_cassettes.push_back(cassettes[i]);
    }
  }
  
  return filtered_cassettes;
}

#endif // COMPLEXITY_FILTER_H