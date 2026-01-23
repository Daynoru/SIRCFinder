// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#ifndef SMART_FILTER_H
#define SMART_FILTER_H

#include "dna_utils.h"
#include <vector>
#include <string>
#include <mutex>
#include <tbb/tbb.h>

using namespace dna_utils;

struct FilterResult {
  bool is_valid_cassette;
  double similarity_score;
};

class CassetteFilter : public RcppParallel::Worker {
private:
  const std::vector<CoreResult>& cassettes_;
  const std::string& genome_sequence_;
  std::vector<FilterResult>& results_;
  double similarity_threshold_;
  mutable tbb::spin_mutex mutex_;
  
public:
  CassetteFilter(const std::vector<CoreResult>& cassettes,
                 const std::string& genome_sequence,
                 std::vector<FilterResult>& results,
                 double similarity_threshold)
    : cassettes_(cassettes),
      genome_sequence_(genome_sequence),
      results_(results),
      similarity_threshold_(similarity_threshold) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      const CoreResult& cassette = cassettes_[i];
      
      if (cassette.positions.empty() || cassette.sequences.empty()) {
        FilterResult result;
        result.is_valid_cassette = true;
        result.similarity_score = 0.0;
        
        tbb::spin_mutex::scoped_lock lock(mutex_);
        results_[i] = result;
        continue;
      }
      
      int dr_length = cassette.sequences[0].length();
      
      std::vector<std::string> spacers = dna_utils::extract_spacers(
        genome_sequence_, cassette.positions, dr_length);
      
      if (spacers.empty()) {
        FilterResult result;
        result.is_valid_cassette = true;
        result.similarity_score = 0.0;
        
        tbb::spin_mutex::scoped_lock lock(mutex_);
        results_[i] = result;
        continue;
      }
      
      bool is_valid = true;
      double max_similarity = 0.0;
      
      for (const std::string& spacer : spacers) {
        if (spacer.empty()) continue;
        
        int lev_dist = weighted_levenshtein(cassette.consensus, spacer);
        
        int window_size = std::min(5, static_cast<int>(spacer.size()));
        int max_matches = sliding_window_match(cassette.consensus, spacer, window_size);
        
        double similarity = calculate_similarity_score(
          lev_dist, max_matches, spacer.size(), window_size);
        
        if (similarity > max_similarity) {
          max_similarity = similarity;
        }
        
        if (similarity > similarity_threshold_) {
          is_valid = false;
          break;
        }
      }
      
      FilterResult result;
      result.is_valid_cassette = is_valid;
      result.similarity_score = max_similarity;
      
      tbb::spin_mutex::scoped_lock lock(mutex_);
      results_[i] = result;
    }
  }
};


inline std::vector<CoreResult> smart_filter_internal(
    const std::vector<CoreResult>& cassettes,
    const std::string& genome_sequence,
    double similarity_threshold,
    int num_threads) {
  
  if (cassettes.empty()) {
    return std::vector<CoreResult>();
  }
  
  std::vector<FilterResult> filter_results(cassettes.size());
  
  CassetteFilter filter(cassettes, genome_sequence, filter_results, similarity_threshold);
  
  try {
    tbb::task_scheduler_init init(num_threads);
    RcppParallel::parallelFor(0, cassettes.size(), filter);
  } catch (const std::exception& e) {
    throw std::runtime_error("Parallel filtering failed: " + std::string(e.what()));
  }
  
  std::vector<CoreResult> filtered_cassettes;
  for (size_t i = 0; i < cassettes.size(); ++i) {
    if (filter_results[i].is_valid_cassette) {
      filtered_cassettes.push_back(cassettes[i]);
    }
  }
  
  return filtered_cassettes;
}

#endif // SMART_FILTER_H