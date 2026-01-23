// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#ifndef EXTENDER_H
#define EXTENDER_H

#include "dna_utils.h"
#include <vector>
#include <string>
#include <tbb/task_scheduler_init.h>

using namespace dna_utils;

// Воркер для обработки ядер
struct CoreExtenderWorker : public RcppParallel::Worker {
  const std::vector<CoreResult>& input_cores;
  const std::string& genome;
  const int max_extension_length;
  const int entropy_window_size;
  std::vector<CoreResult>& output_results;
  
  CoreExtenderWorker(const std::vector<CoreResult>& input_cores,
                     const std::string& genome,
                     std::vector<CoreResult>& output_results,
                     int max_ext, int entropy_win)
    : input_cores(input_cores), 
      genome(genome),
      max_extension_length(max_ext),
      entropy_window_size(entropy_win),
      output_results(output_results) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      output_results[i] = processSingleCore(input_cores[i]);
    }
  }
  
private:
  CoreResult processSingleCore(const CoreResult& core) {
    CoreResult result = core; 
    
    try {
      DrFinderResult right_result = extend_right(core.sequences, core.positions, genome,
                                                 max_extension_length, entropy_window_size);
      
      DrFinderResult final_result = extend_left(right_result.sequences, right_result.positions, genome,
                                                max_extension_length, entropy_window_size);
      
      result.consensus = final_result.consensus;
      result.positions = final_result.positions;
      result.sequences = final_result.sequences;
      result.dr_length = final_result.dr_length; 
      
    } catch (const std::exception& e) {
      
      result.consensus = core.consensus;
    }
    
    return result;
  }
};

inline std::vector<CoreResult> filter_overlapping_cassettes(
    const std::vector<CoreResult>& input_data, int estimated_dr_length) {
  return dna_utils::filter_overlapping_cassettes(input_data, estimated_dr_length);
}

inline std::vector<CoreResult> extend_cores_internal(
    const std::vector<CoreResult>& input_cores,
    const std::string& genome,
    int max_extension_length,
    int entropy_window_size,
    int num_threads) {
  
  if (input_cores.empty()) {
    return std::vector<CoreResult>();
  }
  

  auto filtered_input = input_cores;
  if (filtered_input.empty()) {
    return std::vector<CoreResult>();
  }
  
  std::vector<CoreResult> output_data(filtered_input.size());
  
  try {
    tbb::task_scheduler_init init_process(num_threads);
    CoreExtenderWorker worker(filtered_input, genome, output_data,
                              max_extension_length, entropy_window_size);
    RcppParallel::parallelFor(0, filtered_input.size(), worker);
  } catch (const std::exception& e) {
    throw std::runtime_error("Parallel execution failed: " + std::string(e.what()));
  }
  


  return output_data;
}

#endif // EXTENDER_H