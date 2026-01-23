// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#ifndef SPLIT_CASSETTES_H
#define SPLIT_CASSETTES_H

#include "dna_utils.h"
#include <vector>
#include <string>
#include <algorithm>

using namespace dna_utils;

inline std::vector<CoreResult> split_cassettes_internal(const std::vector<CoreResult>& cassettes) {
  
  if (cassettes.empty()) {
    return std::vector<CoreResult>();
  }
  
  std::vector<CoreResult> result;
  
  for (const auto& cassette : cassettes) {
    if (cassette.positions.size() < 3) {
      continue;
    }
    
    std::vector<int> spacers;
    for (size_t j = 0; j < cassette.positions.size() - 1; j++) {
      int spacer_length = cassette.positions[j + 1] - cassette.positions[j] - cassette.dr_length;
      spacers.push_back(spacer_length);
    }
    
    std::vector<int> split_points;
    int dr_length = cassette.dr_length;
    
    for (size_t j = 0; j < spacers.size(); j++) {
      if (spacers[j] > (2.5 * dr_length) || spacers[j] < (0.5 * dr_length)) {
        split_points.push_back(j); 
      }
    }
    
    if (split_points.empty()) {
      result.push_back(cassette);
      continue;
    }
    
    std::sort(split_points.begin(), split_points.end());
    split_points.push_back(cassette.positions.size() - 1); 
    
    int start_index = 0;
    for (size_t j = 0; j < split_points.size(); j++) {
      int end_index = split_points[j];
      
      if (end_index - start_index + 1 >= 3) {
        CoreResult segment;
        segment.consensus = cassette.consensus;
        segment.dr_length = cassette.dr_length;
        segment.positions.assign(cassette.positions.begin() + start_index, 
                                 cassette.positions.begin() + end_index + 1);
        segment.sequences.assign(cassette.sequences.begin() + start_index, 
                                 cassette.sequences.begin() + end_index + 1);
        
        result.push_back(segment);
      }
      
      start_index = end_index + 1;
    }
  }
  
  return result;
}

#endif // SPLIT_CASSETTES_H