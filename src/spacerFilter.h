// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppParallel)]]

#ifndef SPACER_FILTER_H
#define SPACER_FILTER_H

#include <vector>
#include <string>
#include <algorithm>
#include "dna_utils.h"




 inline bool should_filter_cassette(const std::vector<std::string>& spacers) {
    if (spacers.size() < 2) {
        return false;
    }
    
    int similar_pairs_count = 0;
    int total_pairs = 0;
    
    // Rcpp::Rcout << "=== Checking cassette with " << spacers.size() << " spacers ===" << std::endl;
    
    for (size_t i = 0; i < spacers.size(); ++i) {
        for (size_t j = i + 1; j < spacers.size(); ++j) {
            // Rcpp::Rcout << "Comparing spacer " << i << " vs " << j << std::endl;
            
            bool are_dissimilar = areSequencesDissimilar(spacers[i], spacers[j], 100); 
            
            // Rcpp::Rcout << "Result: " << (are_dissimilar ? "DISSIMILAR" : "SIMILAR") << std::endl;
            
            if (!are_dissimilar) {
                similar_pairs_count++;
            }
            total_pairs++;
        }
    }
    
    // Rcpp::Rcout << "Total pairs: " << total_pairs << ", Similar pairs: " << similar_pairs_count << std::endl;
    // Rcpp::Rcout << "Filter cassette: " << (similar_pairs_count > total_pairs / 2) << std::endl;
    
    return (similar_pairs_count > total_pairs / 2);
}

inline std::vector<CoreResult> spacer_filter_internal(
    const std::vector<CoreResult>& results,
    const std::string& genome_sequence
    ) {
  int total_cassettes = results.size();
  std::vector<CoreResult> filtered_results;
  int filtered_out = 0;
  for (const auto& result : results) {
    std::vector<std::string> spacers = extract_spacers(genome_sequence, result.positions, result.dr_length);
    
    if (!should_filter_cassette(spacers)) {
      filtered_results.push_back(result);
    } else {
      filtered_out++;
    }
  }
  // UNCOMMENT FOR FILTERING STATS
  /*
  Rcpp::Rcout << "==========================================" << std::endl;
  Rcpp::Rcout << "FILTERING STATISTICS:" << std::endl;
  Rcpp::Rcout << "Total cassettes processed: " << total_cassettes << std::endl;
  Rcpp::Rcout << "Cassettes filtered out: " << filtered_out << std::endl;
  Rcpp::Rcout << "Cassettes remaining: " << filtered_results.size() << std::endl;
  Rcpp::Rcout << "Filtering rate: " << (100.0 * filtered_out / total_cassettes) << "%" << std::endl;
  Rcpp::Rcout << "==========================================" << std::endl;
    */

  return filtered_results;
}







#endif // SPACER_FILTER_H