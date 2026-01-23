// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#ifndef CORES_SEARCH_H
#define CORES_SEARCH_H

#include <Rcpp.h>
#include <RcppParallel.h>
#include "dna_utils.h"
#include <RcppParallel/TinyThread.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace dna_utils;

// Parallel cores search struct
struct CoreFinder : public Worker {

  const std::string& genome;
  const int window_size;
  const int step_size;
  const int k_core;
  const int min_matches;
  const int min_repeats;
  const int max_search_offset;
  const int max_extension_length;


  std::vector<std::vector<CoreResult>>& all_results;
  tthread::mutex& mutex;

  CoreFinder(const std::string& genome,
             int window_size,
             int step_size,
             int k_core,
             int min_matches,
             int min_repeats,
             int max_search_offset,
             int max_extension_length,
           std::vector<std::vector<CoreResult>>& all_results,
           tthread::mutex& mutex)
    : genome(genome),
      window_size(window_size),
      step_size(step_size),
      k_core(k_core),
      min_matches(min_matches),
      min_repeats(min_repeats),
      max_search_offset(max_search_offset),
      max_extension_length(max_extension_length),
      all_results(all_results),
      mutex(mutex) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t window_idx = begin; window_idx < end; ++window_idx) {
      process_window(window_idx);
    }
  }

private:
// single window processing
void process_window(std::size_t window_idx) {
    std::vector<CoreResult> window_results;
    std::unordered_set<int> occupied_positions;

    int window_start = window_idx * step_size;
    int window_end = std::min(window_start + window_size, static_cast<int>(genome.size()));
    int first_kmer_end = std::min(window_start + max_search_offset, window_end - k_core);

    for (int pos = window_start; pos < first_kmer_end; ++pos) {
        if (pos + k_core > static_cast<int>(genome.size())) continue;

        // skip occupied positions
        if (occupied_positions.find(pos) != occupied_positions.end()) {
            continue;
        }

        process_position(pos, window_results, occupied_positions, window_start, window_end);
    }




if (!window_results.empty()) {
std::vector<CoreResult> all_split_groups;

    for (const auto& core_group : window_results) {
        auto split_groups = split_core_groups(core_group, k_core, max_extension_length);
        all_split_groups.insert(all_split_groups.end(), split_groups.begin(), split_groups.end());
    }
   // filtering by min_repeats
    std::vector<CoreResult> filtered_results;
    for (const auto& subgroup : all_split_groups) {
        if (subgroup.positions.size() >= static_cast<size_t>(min_repeats)) {
            filtered_results.push_back(subgroup);
        }
    }

    // sorting by positions
    for (size_t i = 0; i < filtered_results.size(); ++i) {
        for (size_t j = i + 1; j < filtered_results.size(); ++j) {
            if (filtered_results[i].positions[0] > filtered_results[j].positions[0]) {
                std::swap(filtered_results[i], filtered_results[j]);
            }
        }
    }

    // saving after filterings
    if (!filtered_results.empty()) {
        save_results(filtered_results);


    }
}

}

// Single genome position
void process_position(int pos,
                     std::vector<CoreResult>& window_results,
                     std::unordered_set<int>& occupied_positions,
                     int window_start,
                     int window_end) {
    std::string core = genome.substr(pos, k_core);

    // Basic filtering
    if (dna_utils::is_softmasked(core)) {
        return;
    }

    if (!dna_utils::is_valid_dna(core)) {
        return;
    }

    if (!dna_utils::crudeFilter(core)) {
        return;
    }

    // All similar cores along max_search_offset
    std::vector<int> matches = find_similar_cores(core, pos);

    // Check min_repeats
    if (matches.size() >= static_cast<size_t>(min_repeats)) {
        CoreResult result;
        result.consensus = core;
        result.dr_length = k_core;
        result.positions = matches;
        result.sequences = extract_sequences(matches);

        window_results.push_back(result);

        // Mark all as occupied positions
        for (int match_pos : matches) {
            occupied_positions.insert(match_pos);

        }
    }
}


  // Similar cores search
  std::vector<int> find_similar_cores(const std::string& core, int start_pos) {
    std::vector<int> matches;
    matches.push_back(start_pos);

    int search_start = start_pos + 1;
    int search_end = std::min(search_start + max_search_offset,
                              static_cast<int>(genome.size()) - k_core);

    for (int search_pos = search_start; search_pos < search_end; ++search_pos) {
      std::string candidate = genome.substr(search_pos, k_core);

      if (is_similar(core, candidate, k_core, min_matches)) {
        matches.push_back(search_pos);
      }
    }

    return matches;
  }

  // Sequences extractions
  std::vector<std::string> extract_sequences(const std::vector<int>& positions) {
    std::vector<std::string> sequences;

    for (int pos : positions) {
      int start = pos;
      if (start + k_core <= static_cast<int>(genome.size())) {
        sequences.push_back(genome.substr(start, k_core));
      }
    }

    return sequences;
  }

  // Saving results
  void save_results(const std::vector<CoreResult>& results) {
    tthread::lock_guard<tthread::mutex> lock(mutex);
    all_results.push_back(results);
  }
};

// Utility for conversion
inline std::vector<CoreResult> flatten_results(const std::vector<std::vector<CoreResult>>& all_results) {
  std::vector<CoreResult> flattened;

  for (const auto& window_results : all_results) {
    for (const auto& core_data : window_results) {
      flattened.push_back(core_data);
    }
  }

  return flattened;
}

// Main cores search function
inline std::vector<CoreResult> deep_cores_search_internal(const std::string& genome,
                                                          int k_core,
                                                          int min_repeats,
                                                          int window_size,
                                                          int step_size,
                                                          int num_threads,
                                                          int min_matches,
                                                          int max_search_offset,
                                                          int max_extension_length,
                                                        bool heuristics) {

  // parameters validation
  if (window_size <= 0 || step_size <= 0 || window_size < k_core) {
    throw std::invalid_argument("Invalid window parameters");
  }

  if (min_matches < 1 || min_matches > k_core) {
    throw std::invalid_argument("min_matches must be between 1 and k_core");
  }

  if (genome.empty()) {
    return std::vector<CoreResult>();
  }

  std::vector<std::vector<CoreResult>> all_results;
  tthread::mutex mutex;

  int num_windows = (genome.length() - window_size + step_size) / step_size;
  if (num_windows <= 0) num_windows = 1;

  CoreFinder finder(genome, window_size, step_size, k_core, min_matches,
                  min_repeats, max_search_offset, max_extension_length,
                  all_results, mutex);

  if (num_threads == 1) {
    finder(0, num_windows);
  } else {
    parallelFor(0, num_windows, finder, num_threads);
  }

  auto flatenned = flatten_results(all_results);

  auto deduplicated = process_final_results(flatenned);

  auto filtered = filter_overlapping_cores_simple(deduplicated, k_core, heuristics);

if (heuristics) {

    for (auto& group : filtered) {
        if (group.dr_length > k_core) {

            std::vector<std::string> extended_sequences;
            for (int pos : group.positions) {

                if (pos + group.dr_length <= genome.length()) {
                    std::string extended_seq = genome.substr(pos, group.dr_length);
                    extended_sequences.push_back(extended_seq);

                } else {

                    extended_sequences.push_back(group.sequences[0]); // fallback
                }
            }


            group.sequences = extended_sequences;
            group.consensus = calculate_consensus(group.sequences);


        }
    }
}

 return filtered;

}

#endif // CORES_SEARCH_H
