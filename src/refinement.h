// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#ifndef REFINEMENT_H
#define REFINEMENT_H

#include "dna_utils.h"
#include <vector>
#include <string>
#include <mutex>
#include <tbb/tbb.h>

using namespace dna_utils;

struct RefinementWorker : public RcppParallel::Worker {
  const std::vector<CoreResult>& input_cassettes;
  const std::string& genome_sequence;
  const double min_similarity;
  const double spacer_similarity_threshold;
  
  std::vector<CoreResult>& output_data;
  
  mutable tbb::spin_mutex output_mutex;
  
  RefinementWorker(const std::vector<CoreResult>& input, 
                   const std::string& genome,
                   double similarity,
                   double spacer_similarity_threshold,
                   std::vector<CoreResult>& output)
    : input_cassettes(input), genome_sequence(genome), 
      min_similarity(similarity), spacer_similarity_threshold(spacer_similarity_threshold), output_data(output) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      CoreResult result = refine_cassette(input_cassettes[i]);
      
      tbb::spin_mutex::scoped_lock lock(output_mutex);
      output_data[i] = result;
    }
  }
  
private:
  std::string extract_sequence(int start, int end) const {
    if (start < 0 || end > static_cast<int>(genome_sequence.size()) || start >= end) {
      return "";
    }
    return genome_sequence.substr(start, end - start);
  }
  
  std::string extract_spacer(int spacer_start, int spacer_end) const {
    return extract_sequence(spacer_start, spacer_end);
  }
  
  bool are_spacers_similar(const std::string& spacer1, const std::string& spacer2, double threshold) const {
    if (spacer1.empty() || spacer2.empty()) return false;
    
    size_t min_len = std::min(spacer1.length(), spacer2.length());
    if (min_len < 5) return false;
    
    int matches = 0;
    for (size_t i = 0; i < min_len; ++i) {
      if (spacer1[i] == spacer2[i]) {
        matches++;
      }
    }
    
    double similarity = static_cast<double>(matches) / min_len;
    return similarity >= threshold;
  }
  
  CoreResult refine_cassette(const CoreResult& cassette) {
    CoreResult result = cassette; 
    bool changed;
    int dr_length = cassette.sequences.empty() ? 0 : cassette.sequences[0].length();
    
    do {
      changed = false;
      
      if (result.positions.empty()) continue;
      
      int first_pos = result.positions.front();
      int search_start_before = std::max(0, first_pos - static_cast<int>(3.5 * dr_length));
      int search_end_before = first_pos - static_cast<int>(1.5 * dr_length);
      
      if (search_end_before > search_start_before) {
        auto before_repeats = find_additional_repeats(genome_sequence, result.consensus,
                                                      dr_length, search_start_before,
                                                      search_end_before, min_similarity);
        
        if (!before_repeats.empty()) {
          std::string new_dr_sequence = extract_sequence(before_repeats.back(), before_repeats.back() + dr_length);
          
          if (!crudeFilter(new_dr_sequence)) {
            continue;
          }
          
          bool is_tandem = false;
          int new_dr_pos = before_repeats.back();
          
          std::string new_spacer = extract_spacer(new_dr_pos + dr_length, first_pos);
          
          if (result.positions.size() >= 2) {
            std::string existing_spacer = extract_spacer(result.positions[0] + dr_length,
                                                         result.positions[1]);
            
            if (are_spacers_similar(new_spacer, existing_spacer, spacer_similarity_threshold)) {
              is_tandem = true;
            }
          }
          
          if (!is_tandem) {
            result.positions.insert(result.positions.begin(), before_repeats.begin(), before_repeats.end());
            
            for (auto it = before_repeats.rbegin(); it != before_repeats.rend(); ++it) {
              std::string seq = extract_sequence(*it, *it + dr_length);
              result.sequences.insert(result.sequences.begin(), seq);
            }
            
            changed = true;
          }
        }
      }
      
      int last_pos = result.positions.back();
      int search_start_after = last_pos + static_cast<int>(1.5 * dr_length);
      int search_end_after = last_pos + static_cast<int>(3.5 * dr_length);
      
      if (search_start_after < static_cast<int>(genome_sequence.length()) - dr_length) {
        auto after_repeats = find_additional_repeats(genome_sequence, result.consensus,
                                                     dr_length, search_start_after,
                                                     search_end_after, min_similarity);
        
        if (!after_repeats.empty()) {
          bool is_tandem = false;
          int new_dr_pos = after_repeats.front();
          
          std::string new_spacer = extract_spacer(last_pos + dr_length, new_dr_pos);
          
          if (result.positions.size() >= 2) {
            int second_last_pos = result.positions[result.positions.size() - 2];
            std::string existing_spacer = extract_spacer(second_last_pos + dr_length, last_pos);
            
            if (are_spacers_similar(new_spacer, existing_spacer, spacer_similarity_threshold)) {
              is_tandem = true;
            }
          }
          
          if (!is_tandem) {
            result.positions.insert(result.positions.end(), after_repeats.begin(), after_repeats.end());
            
            for (const auto& pos : after_repeats) {
              std::string seq = extract_sequence(pos, pos + dr_length);
              result.sequences.push_back(seq);
            }
            
            changed = true;
          }
        }
      }
      
    } while (changed);
    
    return result;
  }
};

inline std::vector<CoreResult> refine_cassettes_internal(
    const std::vector<CoreResult>& cassettes,
    const std::string& genome_sequence,
    double min_similarity,
    double spacer_similarity_threshold,
    int num_threads) {
  
  if (cassettes.empty()) {
    return std::vector<CoreResult>();
  }
  
  std::vector<CoreResult> output_data(cassettes.size());
  
  RefinementWorker worker(cassettes, genome_sequence, min_similarity,spacer_similarity_threshold, output_data);
  
  try {
    tbb::task_scheduler_init init(num_threads);
    
    RcppParallel::parallelFor(0, cassettes.size(), worker);
    
  } catch (const std::exception& e) {
    throw std::runtime_error("Parallel execution failed: " + std::string(e.what()));
  }
  
  return output_data;
}

#endif // REFINEMENT_H