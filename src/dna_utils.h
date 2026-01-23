#ifndef DNA_UTILS_H
#define DNA_UTILS_H

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>
#include <random>
#include <numeric>

namespace dna_utils {

using namespace Rcpp;

inline bool is_softmasked(const std::string& sequence) {
  for (char c : sequence) {
    if (c >= 'a' && c <= 'z') { // Быстрая проверка без islower
      return true;
    }
  }
  return false;
}

inline bool crudeFilter(const std::string& core) {
  if (core.length() < 3) return false;
  
  // Rule 1: minimum 3 unique chars
  std::unordered_set<char> unique_chars;
  for (char c : core) {
    unique_chars.insert(c);
  }
  if (unique_chars.size() < 3) {
    return false;
  }
  
  // Rule 2: no more than 6 ident chars in a row!
  int max_consecutive = 1;
  int current_consecutive = 1;
  char prev_char = core[0];
  
  for (size_t i = 1; i < core.length(); ++i) {
    if (core[i] == prev_char) {
      current_consecutive++;
      max_consecutive = std::max(max_consecutive, current_consecutive);
    } else {
      current_consecutive = 1;
      prev_char = core[i];
    }
    
    if (max_consecutive >= 7) {
      return false;
    }
  }
  
  return true;
}


inline bool is_valid_dna(const std::string& sequence) {
    for (char c : sequence) {
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G' &&
            c != 'a' && c != 't' && c != 'c' && c != 'g') {
            return false;
        }
    }
    return true;
}


inline std::string to_upper_dna(std::string sequence) {
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
    return sequence;
}

inline bool is_similar(const std::string& a, const std::string& b, 
                      int k_core, int min_matches) {
    int mismatches = 0;
    for (int i = 0; i < k_core; ++i) {
        if (a[i] != b[i]) mismatches++;
        if (mismatches > (k_core - min_matches)) return false;
    }
    return true;
}

inline bool already_processed(const std::unordered_set<std::string>& seen_cores, 
                             const std::string& core,
                             int k_core, int min_matches) {
    for (const auto& seen_core : seen_cores) {
        if (is_similar(seen_core, core, k_core, min_matches)) {
            return true;
        }
    }
    return false;
}

struct CoreResult {
    std::string consensus;
    int dr_length;
    std::vector<int> positions;
    std::vector<std::string> sequences;
    void update_sequences(const std::string& genome) {
        sequences.clear();
        for (int pos : positions) {
            if (pos + dr_length <= static_cast<int>(genome.size())) {
                sequences.push_back(genome.substr(pos, dr_length));
            }
        }
    }
};






inline bool has_microsatellite(const std::string& sequence) {
  if (sequence.length() < 6) return false;
  for (size_t i = 0; i <= sequence.length() - 6; ++i) {
    char c = sequence[i];
    if (c == sequence[i+1] && c == sequence[i+2] && 
        c == sequence[i+3] && c == sequence[i+4] && c == sequence[i+5]) {
      return true;
    }
  }
  return false;
}

inline double calculate_position_entropy(const std::vector<std::string>& sequences, size_t pos) {
  if (sequences.empty()) return 0.0;
  
  std::unordered_map<char, int> counts;
  int total = 0;
  
  for (const auto& seq : sequences) {
    if (pos < seq.length()) {
      counts[seq[pos]]++;
      total++;
    }
  }
  
  if (total == 0) return 0.0;
  
  double entropy = 0.0;
  for (const auto& pair : counts) {
    double probability = static_cast<double>(pair.second) / total;
    if (probability > 0) {
      entropy -= probability * log2(probability);
    }
  }
  
  return entropy;
}

inline double calculate_region_entropy(const std::vector<std::string>& sequences, 
                                       size_t start_pos, size_t end_pos) {
  if (sequences.empty() || start_pos >= end_pos) return 0.0;
  
  double total_entropy = 0.0;
  size_t count = 0;
  
  for (size_t pos = start_pos; pos < end_pos; ++pos) {
    double pos_entropy = calculate_position_entropy(sequences, pos);
    if (pos_entropy > 0) {
      total_entropy += pos_entropy;
      count++;
    }
  }
  
  return count > 0 ? total_entropy / count : 0.0;
}

struct EntropyResult {
  double dr_entropy;
  double spacer_entropy;
  double score;
  int extension_length;
  
  EntropyResult() : dr_entropy(0.0), spacer_entropy(0.0), 
  score(-std::numeric_limits<double>::infinity()), extension_length(0) {}
};


inline List convert_results_to_list(const std::vector<std::vector<CoreResult>>& all_results,
                                    const std::string& genome, int k_core) {
  size_t total_cores = 0;
  for (const auto& window_results : all_results) {
    total_cores += window_results.size();
  }
  
  List results(total_cores);
  size_t index = 0;
  
  for (const auto& window_results : all_results) {
    for (const auto& core_data : window_results) {
      List core_list;
      core_list["consensus"] = core_data.consensus;
      core_list["positions"] = core_data.positions;
      
      std::vector<std::string> actual_sequences;
      for (int pos : core_data.positions) {
        int start = pos - 1; 
        if (start + k_core <= static_cast<int>(genome.size())) {
          actual_sequences.push_back(genome.substr(start, k_core));
        }
      }
      
      core_list["sequences"] = actual_sequences;
      results[index++] = core_list;
    }
  }
  
  return results;
}



double calcEntropy(const std::vector<std::string>& sequences) {
  if (sequences.empty()) return 0.0;
  
  int seq_length = sequences[0].length();
  double total_entropy = 0.0;
  int valid_columns = 0;
  
  for (int col = 0; col < seq_length; ++col) {
    std::map<char, int> nucleotide_counts = {{'A', 0}, {'T', 0}, {'C', 0}, {'G', 0}};
    int total_nucleotides = 0;
    
    for (size_t row = 0; row < sequences.size(); ++row) {
      if (col >= sequences[row].length()) continue;
      
      char nucleotide = toupper(sequences[row][col]);
      if (nucleotide == 'A' || nucleotide == 'T' || nucleotide == 'C' || nucleotide == 'G') {
        nucleotide_counts[nucleotide]++;
        total_nucleotides++;
      }
    }
    
    if (total_nucleotides == 0) continue;
    
    double col_entropy = 0.0;
    for (auto const& pair : nucleotide_counts) {
      if (pair.second > 0) {
        double p = (double)pair.second / total_nucleotides;
        col_entropy -= p * log2(p);
      }
    }
    
    total_entropy += col_entropy;
    valid_columns++;
  }
  
  return valid_columns == 0 ? 0.0 : total_entropy / valid_columns;
}

struct DrFinderResult {
  std::string consensus;
  int dr_length;
  std::vector<int> positions; // 0-based 
  std::vector<std::string> sequences;
};

std::string calculate_consensus(const std::vector<std::string>& sequences) {
  if (sequences.empty()) return "";
  
  size_t length = sequences[0].length();
  std::string consensus;
  consensus.reserve(length);
  
  for (size_t i = 0; i < length; ++i) {
    int countA = 0, countC = 0, countG = 0, countT = 0;
    
    for (const auto& seq : sequences) {
      if (i >= seq.length()) continue;
      
      switch (seq[i]) {
      case 'A': countA++; break;
      case 'C': countC++; break;
      case 'G': countG++; break;
      case 'T': countT++; break;
      }
    }
    
    int max_count = std::max({countA, countC, countG, countT});
    if (max_count == countA) consensus += 'A';
    else if (max_count == countC) consensus += 'C';
    else if (max_count == countG) consensus += 'G';
    else consensus += 'T';
  }
  
  return consensus;
}


inline std::pair<double, int> find_boundary_by_entropy_jump(const std::vector<std::string>& alignment, 
                                    int boundary_position) {
    std::vector<double> positional_entropy;
    for (size_t pos = 0; pos < alignment[0].length(); ++pos) {
        positional_entropy.push_back(calculate_position_entropy(alignment, pos));
    }
    
    int num_sequences = alignment.size();
    int k = (num_sequences + 1) / 2;
    int base_mutations = k - 1;
    
    double base_threshold = 0.0;
    if (base_mutations > 0 && base_mutations < num_sequences) {
        double p_major = (num_sequences - base_mutations) / (double)num_sequences;
        double p_minor = 1.0 / num_sequences;
        base_threshold = -p_major * log2(p_major) - base_mutations * p_minor * log2(p_minor);
    }
    
    double max_repeat_entropy = 0.0;
    int repeat_region_length = std::min(10, static_cast<int>(positional_entropy.size()));
    
    for (int pos = 0; pos < repeat_region_length; ++pos) {
        if (positional_entropy[pos] > max_repeat_entropy) {
            max_repeat_entropy = positional_entropy[pos];
        }
    }
    
    double final_threshold = base_threshold;
    bool used_adaptive = false;
    if (max_repeat_entropy > base_threshold * 0.8) {
        used_adaptive = true;
        if (k > 0 && k < num_sequences) {
            double p_major = (num_sequences - k) / (double)num_sequences;
            double p_minor = 1.0 / num_sequences;
            final_threshold = -p_major * log2(p_major) - k * p_minor * log2(p_minor);
        }
    }
    
    int found_boundary = -1;
    for (size_t pos = repeat_region_length; pos < positional_entropy.size(); ++pos) {
        if (positional_entropy[pos] > final_threshold) {
            found_boundary = pos;
            break;
        }
    }
    
   
    
    if (found_boundary != -1) {
        double score = positional_entropy[found_boundary] - max_repeat_entropy;
        return std::make_pair(score, found_boundary);
    }
    
    return std::make_pair(0.0, -1);
}

DrFinderResult extend_right(const std::vector<std::string>& input_sequences,
                            const std::vector<int>& input_positions,
                            const std::string& genome,
                            int max_extension_length,
                            int original_dr_length) {
  
  DrFinderResult best_result;
  best_result.consensus = "";
  best_result.dr_length = original_dr_length;
  best_result.positions = input_positions;
  best_result.sequences = input_sequences;
  
  std::vector<int> base_positions = input_positions;
  int original_length = input_sequences.empty() ? 0 : input_sequences[0].length();
  
  int min_spacer_length = std::numeric_limits<int>::max();
  for (size_t i = 0; i < base_positions.size() - 1; ++i) {
    int spacer_length = base_positions[i + 1] - base_positions[i] - original_length;
    if (spacer_length < min_spacer_length) {
      min_spacer_length = spacer_length;
    }
  }
  
  int last_repeat_end = base_positions.back() + original_length;
  int genomic_remaining = genome.length() - last_repeat_end;
  
  int max_possible_extension = std::min({min_spacer_length, genomic_remaining, max_extension_length});
  
  if (max_possible_extension > 0) {
    int analysis_length = 10 + max_possible_extension; 
    
    std::vector<std::string> full_alignment;
    for (size_t i = 0; i < base_positions.size(); ++i) {
      int repeat_end = base_positions[i] + original_length;
      int analysis_start = repeat_end - 10;
      full_alignment.push_back(genome.substr(analysis_start, analysis_length));
    }
    
    auto result = find_boundary_by_entropy_jump(full_alignment, 10);
    int boundary_position = result.second;
    
    int optimal_extension = 0;
    if (boundary_position != -1) {
      optimal_extension = boundary_position - 10;
      optimal_extension = std::min(optimal_extension, max_possible_extension);
      optimal_extension = std::max(optimal_extension, 0);
    }
    
    int best_length = original_length + optimal_extension;
    if (optimal_extension >= 0) {
      std::vector<std::string> repeats_only;
      bool all_valid = true;
      
      for (size_t i = 0; i < base_positions.size(); ++i) {
        int start = base_positions[i];
        int end = start + best_length;
        if (start >= 0 && end <= static_cast<int>(genome.length())) {
          repeats_only.push_back(genome.substr(start, best_length));
        } else {
          all_valid = false;
          break;
        }
      }
      
      if (all_valid && !repeats_only.empty()) {
        best_result.consensus = calculate_consensus(repeats_only);
        best_result.dr_length = best_length;
        best_result.sequences = repeats_only;
      }
    }
  }
  
  return best_result;
}




inline std::pair<double, int> find_boundary_by_entropy_jump_LEFT(const std::vector<std::string>& alignment, 
                                    int boundary_position) {
    std::vector<double> positional_entropy;
    for (size_t pos = 0; pos < alignment[0].length(); ++pos) {
        positional_entropy.push_back(calculate_position_entropy(alignment, pos));
    }
    
    int num_sequences = alignment.size();
    int k = (num_sequences + 1) / 2;
    int base_mutations = std::max(1, k - 2);
    
    double base_threshold = 0.0;
    if (base_mutations > 0 && base_mutations < num_sequences) {
        double p_major = (num_sequences - base_mutations) / (double)num_sequences;
        double p_minor = 1.0 / num_sequences;
        base_threshold = -p_major * log2(p_major) - base_mutations * p_minor * log2(p_minor);
    }
    
    double max_repeat_entropy = 0.0;
    int alignment_length = alignment[0].length();
    int repeat_region_start = std::max(0, alignment_length - 10); // последние 10 позиций
    int repeat_region_length = alignment_length - repeat_region_start;
    
    for (int pos = repeat_region_start; pos < alignment_length; ++pos) {
        if (positional_entropy[pos] > max_repeat_entropy) {
            max_repeat_entropy = positional_entropy[pos];
        }
    }
    
    double final_threshold = base_threshold;
    bool used_adaptive = false;


    
    int found_boundary = -1;
    for (int pos = repeat_region_start - 1; pos >= 0; --pos) {
        if (positional_entropy[pos] >= final_threshold) {
            found_boundary = pos;
            break;
        }
    }
    if (found_boundary != -1) {
        double score = positional_entropy[found_boundary] - max_repeat_entropy;
        return std::make_pair(score, found_boundary);
    }
    
    return std::make_pair(0.0, -1);
}







DrFinderResult extend_left(const std::vector<std::string>& input_sequences,
                           const std::vector<int>& input_positions,
                           const std::string& genome,
                           int max_extension_length,
                           int original_dr_length) {
  
  DrFinderResult best_result;
  best_result.consensus = "";
  best_result.dr_length = original_dr_length;
  best_result.positions = input_positions;
  best_result.sequences = input_sequences;
  
  std::vector<int> base_positions = input_positions;
  int original_length = input_sequences.empty() ? 0 : input_sequences[0].length();
  
  int min_spacer_length = std::numeric_limits<int>::max();
  for (size_t i = 0; i < base_positions.size() - 1; ++i) {
    int spacer_length = base_positions[i + 1] - base_positions[i] - original_length;
    if (spacer_length < min_spacer_length) {
      min_spacer_length = spacer_length;
    }
  }
  
  int first_repeat_start = base_positions.front();
  int genomic_remaining_left = first_repeat_start; 
  
  int max_possible_extension = std::min({min_spacer_length, genomic_remaining_left, max_extension_length});
  
  if (max_possible_extension > 0) {
    int analysis_length = 10 + max_possible_extension; // 10 нт повтора + спейсер слева
    
    std::vector<std::string> full_alignment;
    for (size_t i = 0; i < base_positions.size(); ++i) {
      int repeat_start = base_positions[i];
      int analysis_start = repeat_start - max_possible_extension;
      full_alignment.push_back(genome.substr(analysis_start, analysis_length));
    }
    
    auto result = find_boundary_by_entropy_jump_LEFT(full_alignment, max_possible_extension);
    int boundary_position = result.second;
    
    int optimal_extension = 0;
    if (boundary_position != -1) {
      optimal_extension =  boundary_position - max_possible_extension;
      optimal_extension = std::min(optimal_extension, max_possible_extension);
      optimal_extension = std::max(optimal_extension, 0);
    }
    
    int best_length = original_length + optimal_extension;
    if (optimal_extension >= 0) {
      std::vector<std::string> repeats_only;
      std::vector<int> new_positions;
      bool all_valid = true;
      
      for (size_t i = 0; i < base_positions.size(); ++i) {
        int start = base_positions[i] - optimal_extension;
        int end = start + best_length;
        if (start >= 0 && end <= static_cast<int>(genome.length())) {
          repeats_only.push_back(genome.substr(start, best_length));
          new_positions.push_back(start);
        } else {
          all_valid = false;
          break;
        }
      }
      
      if (all_valid && !repeats_only.empty()) {
        best_result.consensus = calculate_consensus(repeats_only);
        best_result.dr_length = best_length;
        best_result.positions = new_positions; // ОБНОВЛЯЕМ ПОЗИЦИИ!
        best_result.sequences = repeats_only;
      }
    }
  }
  
  return best_result;
}


template<typename T>
std::vector<T> filter_overlapping_cassettes(const std::vector<T>& input_data, int dr_length) {
  std::vector<T> filtered_data;
  if (input_data.empty()) return filtered_data;

  std::vector<bool> keep(input_data.size(), true);

  for (size_t i = 0; i < input_data.size(); ++i) {
    if (!keep[i]) continue;

    for (size_t j = i + 1; j < input_data.size(); ++j) {
      if (!keep[j]) continue;

      bool overlaps = false;
      for (size_t k = 0; k < input_data[i].positions.size() && !overlaps; ++k) {
        for (size_t l = 0; l < input_data[j].positions.size(); ++l) {
          if (std::abs(input_data[i].positions[k] - input_data[j].positions[l]) < dr_length) {
            overlaps = true;
            break;
          }
        }
      }

      if (overlaps) {
        if (input_data[i].positions.size() >= input_data[j].positions.size()) {
          keep[j] = false;
        } else {
          keep[i] = false;
          break;
        }
      }
    }
  }

  for (size_t i = 0; i < input_data.size(); ++i) {
    if (keep[i]) {
      filtered_data.push_back(input_data[i]);
    }
  }

  return filtered_data;
}

struct Cassette {
  std::string consensus;
  int dr_length;
  std::vector<int> positions;
  std::vector<std::string> sequences;
};

inline std::vector<int> calculate_spacers(const std::vector<int>& positions, int dr_length) {
  std::vector<int> spacers;
  if (positions.size() < 2) {
    return spacers;
  }
  
  for (size_t i = 0; i < positions.size() - 1; i++) {
    int spacer_length = positions[i + 1] - positions[i] - dr_length;
    spacers.push_back(spacer_length);
  }
  
  return spacers;
}

inline bool should_split_spacer(int spacer_length, int dr_length) {
  return spacer_length < (0.5 * dr_length) || spacer_length > (2.5 * dr_length);
}

inline std::vector<Cassette> split_cassette(const Cassette& original, size_t split_index) {
  std::vector<Cassette> result;
  
  Cassette first_part;
  first_part.consensus = original.consensus;
  first_part.dr_length = original.dr_length;
  
  for (size_t i = 0; i <= split_index; i++) {
    first_part.positions.push_back(original.positions[i]);
    first_part.sequences.push_back(original.sequences[i]);
  }
  
  Cassette second_part;
  second_part.consensus = original.consensus;
  second_part.dr_length = original.dr_length;
  
  for (size_t i = split_index + 1; i < original.positions.size(); i++) {
    second_part.positions.push_back(original.positions[i]);
    second_part.sequences.push_back(original.sequences[i]);
  }
  
  result.push_back(first_part);
  result.push_back(second_part);
  
  return result;
}



inline double sequence_similarity(const std::string& seq, const std::string& consensus) {
  if (seq.length() != consensus.length() || seq.empty()) {
    return 0.0;
  }
  
  int matches = 0;
  for (size_t i = 0; i < seq.length(); i++) {
    if (seq[i] == consensus[i]) {
      matches++;
    }
  }
  
  return static_cast<double>(matches) / seq.length();
}

inline std::string extract_sequence(const std::string& genome, int start, int end) {
  if (start < 0) start = 0;
  if (end > static_cast<int>(genome.length())) end = genome.length();
  if (start >= end) return "";
  
  return genome.substr(start, end - start);
}

inline bool is_valid_repeat(const std::string& sequence, const std::string& consensus, 
                            double min_similarity = 0.8) {
  if (sequence.length() != consensus.length()) return false;
  double similarity = sequence_similarity(sequence, consensus);
  return similarity >= min_similarity;
}

inline std::vector<int> find_additional_repeats(const std::string& genome, 
                                                const std::string& consensus,
                                                int dr_length,
                                                int search_start, 
                                                int search_end,
                                                double min_similarity = 0.8) {
  std::vector<int> found_positions;
  
  for (int pos = search_start; pos <= search_end - dr_length; pos++) {
    std::string candidate = extract_sequence(genome, pos, pos + dr_length);
    if (is_valid_repeat(candidate, consensus, min_similarity)) {
      found_positions.push_back(pos);
    }
  }
  
  return found_positions;
}


inline int weighted_levenshtein(const std::string& s1, const std::string& s2, 
                                int sub_cost = 1, int indel_cost = 3) {
  size_t len1 = s1.size();
  size_t len2 = s2.size();
  
  std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1));
  
  for (size_t i = 0; i <= len1; ++i) {
    dp[i][0] = i * indel_cost;
  }
  for (size_t j = 0; j <= len2; ++j) {
    dp[0][j] = j * indel_cost;
  }
  
  for (size_t i = 1; i <= len1; ++i) {
    for (size_t j = 1; j <= len2; ++j) {
      if (s1[i-1] == s2[j-1]) {
        dp[i][j] = dp[i-1][j-1];
      } else {
        int substitution = dp[i-1][j-1] + sub_cost;
        int deletion = dp[i-1][j] + indel_cost;
        int insertion = dp[i][j-1] + indel_cost;
        dp[i][j] = std::min({substitution, deletion, insertion});
      }
    }
  }
  
  return dp[len1][len2];
}

inline int sliding_window_match(const std::string& consensus, 
                                const std::string& spacer, 
                                int window_size) {
  int max_matches = 0;
  int spacer_len = spacer.size();
  int consensus_len = consensus.size();
  
  if (spacer_len < window_size || consensus_len < window_size) {
    return 0;
  }
  
  for (int i = 0; i <= spacer_len - window_size; ++i) {
    std::string spacer_window = spacer.substr(i, window_size);
    
    for (int j = 0; j <= consensus_len - window_size; ++j) {
      std::string consensus_window = consensus.substr(j, window_size);
      
      int matches = 0;
      for (int k = 0; k < window_size; ++k) {
        if (spacer_window[k] == consensus_window[k]) {
          matches++;
        }
      }
      
      if (matches > max_matches) {
        max_matches = matches;
      }
    }
  }
  
  return max_matches;
}

inline std::vector<std::string> extract_spacers(const std::string& genome, 
                                                const std::vector<int>& positions, 
                                                int dr_length) {
  std::vector<std::string> spacers;
  
  if (positions.size() < 2) {
    return spacers;
  }
  
  for (size_t i = 0; i < positions.size() - 1; ++i) {
    int spacer_start = positions[i] + dr_length;
    int spacer_end = positions[i + 1];
    
    if (spacer_start >= 0 && spacer_end > spacer_start && spacer_end <= genome.size()) {
      int length = spacer_end - spacer_start;
      if (length > 0) {
        std::string spacer = genome.substr(spacer_start, length);
        spacers.push_back(spacer);
      }
    }
  }
  
  return spacers;
}

inline double calculate_similarity_score(int levenshtein_dist, int max_matches, 
                                         int spacer_length, int window_size = 5) {
  double score1 = 1.0 - (static_cast<double>(levenshtein_dist) / 
                         (spacer_length * 3.0)); 
  
  double score2 = static_cast<double>(max_matches) / window_size; 
  
  return (score1 + score2) / 2.0; 
}


inline std::string extract_spacer(const std::string& genome, int dr1_end, int dr2_start) {
  if (dr1_end >= dr2_start) {
    return "";
  }
  return genome.substr(dr1_end, dr2_start - dr1_end);
}

inline bool are_spacers_similar(const std::string& spacer1, const std::string& spacer2, 
                                double similarity_threshold = 0.8) {
  if (spacer1.empty() || spacer2.empty()) {
    return false;
  }
  
  int max_len = std::max(spacer1.length(), spacer2.length());
  if (max_len == 0) {
    return true;
  }
  
  int distance = weighted_levenshtein(spacer1, spacer2, 1, 3);
  double similarity = 1.0 - (static_cast<double>(distance) / max_len);
  
  return similarity >= similarity_threshold;
}


inline double calculate_complexity(const std::string& sequence, int k) {
  if (sequence.length() < k) return 0.0;
  
  std::unordered_set<std::string> kmers;
  for (size_t i = 0; i <= sequence.length() - k; ++i) {
    kmers.insert(sequence.substr(i, k));
  }
  
  int unique_kmers = kmers.size();
  int max_possible_kmers = sequence.length() - k + 1;
  
  return static_cast<double>(unique_kmers) / max_possible_kmers;
}

inline bool passes_complexity_threshold(const std::string& sequence, int max_k, double threshold) {
  if (sequence.empty()) return false;
  
  for (int k = 3; k <= max_k; ++k) {
    if (sequence.length() >= k) {
      double complexity = calculate_complexity(sequence, k);
      if (complexity < threshold) {
        return false;
      }
    }
  }
  return true;
}



std::vector<CoreResult> split_core_groups(const CoreResult& core_data, 
                                         int k_core, 
                                         int max_extension_length) {
    std::vector<CoreResult> sub_groups;
    
    if (core_data.positions.empty()) return sub_groups;
    
    std::vector<int> sorted_positions = core_data.positions;
    std::vector<std::string> sorted_sequences = core_data.sequences;
    
    for (size_t i = 0; i < sorted_positions.size(); ++i) {
        for (size_t j = i + 1; j < sorted_positions.size(); ++j) {
            if (sorted_positions[i] > sorted_positions[j]) {
                std::swap(sorted_positions[i], sorted_positions[j]);
                std::swap(sorted_sequences[i], sorted_sequences[j]);
            }
        }
    }
    
    int min_distance = k_core + static_cast<int>(0.5 * k_core); // 15
    int max_distance = max_extension_length + static_cast<int>(2.5 * max_extension_length); // 193
    
    std::vector<int> current_subgroup_positions;
    std::vector<std::string> current_subgroup_sequences;
    
    current_subgroup_positions.push_back(sorted_positions[0]);
    current_subgroup_sequences.push_back(sorted_sequences[0]);
    
    for (size_t i = 1; i < sorted_positions.size(); ++i) {
        int last_pos = current_subgroup_positions.back();
        int current_pos = sorted_positions[i];
        int actual_distance = current_pos - last_pos;
        
        if (actual_distance >= min_distance && actual_distance <= max_distance) {
            current_subgroup_positions.push_back(current_pos);
            current_subgroup_sequences.push_back(sorted_sequences[i]);
        } else {
            if (current_subgroup_positions.size() >= 1) {
                CoreResult new_group = core_data;
                new_group.positions = current_subgroup_positions;
                new_group.sequences = current_subgroup_sequences;
                sub_groups.push_back(new_group);
            }
            current_subgroup_positions.clear();
            current_subgroup_sequences.clear();
            current_subgroup_positions.push_back(current_pos);
            current_subgroup_sequences.push_back(sorted_sequences[i]);
        }
    }
    
    if (!current_subgroup_positions.empty()) {
        CoreResult new_group = core_data;
        new_group.positions = current_subgroup_positions;
        new_group.sequences = current_subgroup_sequences;
        new_group.consensus = calculate_consensus(current_subgroup_sequences);
        sub_groups.push_back(new_group);
    }
    
    return sub_groups;
}

inline double calculate_total_entropy(const std::vector<std::string>& sequences, 
                                      size_t start_pos, size_t end_pos) {
  if (sequences.empty() || start_pos >= end_pos) return 0.0;
  
  double total_entropy = 0.0;
  
  for (size_t pos = start_pos; pos < end_pos; ++pos) {
    total_entropy += calculate_position_entropy(sequences, pos);
  }
  
  return total_entropy;
}


inline bool are_cores_identical(const CoreResult& a, const CoreResult& b) {
    return a.consensus == b.consensus && 
           a.positions == b.positions && 
           a.sequences == b.sequences;
}

inline std::vector<CoreResult> process_final_results(std::vector<CoreResult> flattened) {
    std::sort(flattened.begin(), flattened.end(), 
        [](const CoreResult& a, const CoreResult& b) {
            return a.positions[0] < b.positions[0];
        });
    
    auto last = std::unique(flattened.begin(), flattened.end(), are_cores_identical);
    flattened.erase(last, flattened.end());
    
    return flattened;
}

inline bool do_groups_overlap(const CoreResult& a, const CoreResult& b, int k_core) {
    for (int pos_a : a.positions) {
        for (int pos_b : b.positions) {
            if (std::abs(pos_a - pos_b) < k_core) {
                return true; 
            }
        }
    }
    return false;
}




template<typename T>
inline std::vector<T> select_best_from_clusters(const std::vector<T>& input_data,
                                        const std::vector<double>& weights,
                                        const std::vector<int>& overlap_clusters) {
    int max_cluster = *std::max_element(overlap_clusters.begin(), overlap_clusters.end());
    
    std::vector<T> result;
    
    for (int cluster_id = 0; cluster_id <= max_cluster; ++cluster_id) {
        int best_index = -1;
        double best_weight = -1.0;
        size_t best_repeats = 0;
        
        for (size_t i = 0; i < input_data.size(); ++i) {
            if (overlap_clusters[i] == cluster_id) {


                if (weights[i] > best_weight || 
                   (weights[i] == best_weight && input_data[i].positions.size() > best_repeats)) {
                    best_weight = weights[i];
                    best_repeats = input_data[i].positions.size();
                    best_index = i;
                }
            }
        }
        
        if (best_index != -1) {
            result.push_back(input_data[best_index]);
        }
    }
    
    return result;
}



template<typename T>
inline std::vector<T> heuristic_selection(const std::vector<T>& input_data, 
                                         const std::vector<double>& weights,
                                         const std::vector<int>& overlap_clusters,
                                         int k_core) {
    std::vector<T> result;
    std::vector<bool> processed(input_data.size(), false);
    
    for (size_t i = 0; i < input_data.size(); ++i) {
        if (processed[i]) continue;
        
        int cluster_id = overlap_clusters[i];
        
        std::vector<size_t> cluster_indices;
        for (size_t j = 0; j < input_data.size(); ++j) {
            if (overlap_clusters[j] == cluster_id && !processed[j]) {
                cluster_indices.push_back(j);
            }
        }
        
        if (cluster_indices.empty()) continue;
        
        auto chain = find_heuristic_chain(input_data, weights, cluster_indices, k_core);
        
        if (chain.size() > 1) {
            T merged = merge_heuristic_chain(input_data, chain, k_core);
            result.push_back(merged);
        } else {
            size_t best_index = cluster_indices[0];
            double best_weight = weights[best_index];
            
            for (size_t idx : cluster_indices) {
                if (weights[idx] > best_weight) {
                    best_weight = weights[idx];
                    best_index = idx;
                }
            }
            result.push_back(input_data[best_index]);
        }
        
        for (size_t idx : cluster_indices) {
            processed[idx] = true;
        }
    }
    
    return result;
}

template<typename T>
inline std::vector<size_t> find_heuristic_chain(const std::vector<T>& input_data,
                                               const std::vector<double>& weights,
                                               const std::vector<size_t>& cluster_indices,
                                               int k_core) {
    std::vector<size_t> chain;
    
    std::vector<size_t> perfect_groups;
    for (size_t idx : cluster_indices) {
        if (weights[idx] >= 0.999) {
            perfect_groups.push_back(idx);
        }
    }
    
    if (perfect_groups.size() < 2) return chain;
    
    std::sort(perfect_groups.begin(), perfect_groups.end(),
              [&](size_t a, size_t b) {
                  return input_data[a].positions[0] < input_data[b].positions[0];
              });
    
    chain.push_back(perfect_groups[0]);
    const auto& first_group = input_data[perfect_groups[0]];
    
    for (size_t i = 1; i < perfect_groups.size(); ++i) {
        const auto& current_group = input_data[perfect_groups[i]];
        const auto& prev_group = input_data[chain.back()];
        
        int shift = current_group.positions[0] - prev_group.positions[0];
        if (shift >= 2*k_core + 1 || shift <= 0) break;
        
        bool constant_shift = true;
        if (current_group.positions.size() == prev_group.positions.size()) {
            for (size_t j = 0; j < current_group.positions.size(); ++j) {
                if (current_group.positions[j] - prev_group.positions[j] != shift) {
                    constant_shift = false;
                    break;
                }
            }
        } else {
            constant_shift = false;
        }
        
        if (constant_shift) {
            chain.push_back(perfect_groups[i]);
        } else {
            break;
        }
    }
    
    return chain;
}

template<typename T>
inline T merge_heuristic_chain(const std::vector<T>& input_data,
                              const std::vector<size_t>& chain,
                              int k_core) {
    T merged = input_data[chain[0]];
    
    if (chain.size() > 1) {
        int last_pos = input_data[chain.back()].positions[0];
        merged.dr_length = last_pos + k_core - merged.positions[0];
        
       
    }
    
    return merged;
}


template<typename T>
inline std::vector<T> filter_overlapping_cores_simple(const std::vector<T>& input_data, int k_core, bool heuristics) {
    if (input_data.empty()) return input_data;

  
    std::vector<double> weights;
    for (const auto& core : input_data) {
        double Es = calculate_total_entropy(core.sequences, 0, k_core);
        double weight = (2.0 * k_core - Es) / (2.0 * k_core);
        weights.push_back(weight);
    }
    
    std::vector<int> overlap_clusters = cluster_overlapping_groups_simple(input_data, k_core);
    
    if (!heuristics) {
    return select_best_from_clusters(input_data, weights, overlap_clusters);
} else {
    auto heuristic_results = heuristic_selection(input_data, weights, overlap_clusters, k_core);
    
   
    
    return heuristic_results;
}
}

struct Interval {
    int start;
    int end;
    int original_index;
};

inline std::vector<int> cluster_overlapping_groups_simple(const std::vector<CoreResult>& data, int k_core) {
    if (data.empty()) return std::vector<int>();
    
    std::vector<Interval> intervals;
    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i].positions.empty()) continue;
        int start = data[i].positions[0];
        int end = data[i].positions.back() + k_core;
        intervals.push_back({start, end, static_cast<int>(i)});
    }
    
    std::vector<Interval> reduced;
    if (!intervals.empty()) {
        Interval current = intervals[0];
        
        for (size_t i = 1; i < intervals.size(); ++i) {
            if (intervals[i].start <= current.end) {
                current.end = std::max(current.end, intervals[i].end);
            } else {
                reduced.push_back(current);
                current = intervals[i];
            }
        }
        reduced.push_back(current);
    }

    
    std::vector<int> clusters(data.size(), -1);
    
    for (size_t i = 0; i < data.size(); ++i) {
        int start = data[i].positions[0];
        int end = data[i].positions.back() + k_core;
        
        for (size_t j = 0; j < reduced.size(); ++j) {
            if (start <= reduced[j].end && end >= reduced[j].start) {
                clusters[i] = j;
                break;
            }
        }
    }
    
    return clusters;
}

//Needleman-Wunsch

const int MATCH_SCORE = 5;
const int MISMATCH_SCORE = -4; 
const int GAP_PENALTY = -8;

double calculatePID(const std::string& seq1, const std::string& seq2);
std::string shuffleSequence(const std::string& sequence);
bool areSequencesDissimilar(const std::string& seq1, const std::string& seq2, int num_shuffles);




inline double needleman_wunsch_identity(const std::string& seq1, const std::string& seq2, 
                                int gap_penalty = -1, int mismatch_score = -1, int match_score = 1) {
    int n = seq1.length();
    int m = seq2.length();
    
    std::vector<std::vector<int>> matrix(n + 1, std::vector<int>(m + 1, 0));
    
    for (int i = 0; i <= n; i++) 
        matrix[i][0] = i * gap_penalty;
    for (int j = 0; j <= m; j++) 
        matrix[0][j] = j * gap_penalty;
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int score = (seq1[i-1] == seq2[j-1]) ? match_score : mismatch_score;
            int match = matrix[i-1][j-1] + score;
            int del = matrix[i-1][j] + gap_penalty;
            int ins = matrix[i][j-1] + gap_penalty;
            matrix[i][j] = std::max({match, del, ins});
        }
    }
    
    std::string aligned1, aligned2;
    int i = n, j = m;
    int max_steps = n + m + 100;
    
    while ((i > 0 || j > 0) && max_steps-- > 0) {
        int score = (i > 0 && j > 0) ? ((seq1[i-1] == seq2[j-1]) ? match_score : mismatch_score) : 0;
        
        if (i > 0 && j > 0 && matrix[i][j] == matrix[i-1][j-1] + score) {
            aligned1 = seq1[i-1] + aligned1;
            aligned2 = seq2[j-1] + aligned2;
            i--; j--;
        }
        else if (i > 0 && matrix[i][j] == matrix[i-1][j] + gap_penalty) {
            aligned1 = seq1[i-1] + aligned1;
            aligned2 = '-' + aligned2;
            i--;
        }
        else if (j > 0) {
            aligned1 = '-' + aligned1;
            aligned2 = seq2[j-1] + aligned2;
            j--;
        }
        else {
            break;
        }
    }
    
    if (aligned1.empty()) return 0.0;
    
    int matches = 0, compared = 0;
    for (size_t i = 0; i < aligned1.length(); i++) {
        if (aligned1[i] != '-' && aligned2[i] != '-') {
            compared++;
            if (aligned1[i] == aligned2[i]) matches++;
        }
    }
    
    return (compared > 0) ? (matches * 100.0 / compared) : 0.0;
}





std::string shuffleSequence(const std::string& sequence) {
    int length = sequence.size();  
    std::string bases = "ACGT";
    std::string result;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 3);
    
    for (int i = 0; i < length; i++) {
        result += bases[dis(gen)];
    }
    return result;
}

bool areSequencesDissimilar(const std::string& seq1, const std::string& seq2, int num_shuffles = 1000) {
    
  const double SIMILARITY_THRESHOLD = 90.0;
    double original_pid = needleman_wunsch_identity(seq1, seq2);

    if (original_pid > SIMILARITY_THRESHOLD) {
                return false;
            }
    
    std::vector<double> shuffle_pids;
    for (int i = 0; i < num_shuffles; i++) {
        std::string shuffled1 = shuffleSequence(seq1);
        std::string shuffled2 = shuffleSequence(seq2);
        double shuffle_pid = needleman_wunsch_identity(shuffled1, shuffled2);
        shuffle_pids.push_back(shuffle_pid);
    }
    
    // 3. Расчет p-value
    int count_better = 0;
    for (double shuffle_pid : shuffle_pids) {
        if (shuffle_pid >= original_pid) {
            count_better++;
        }
    }
    double p_value = static_cast<double>(count_better) / num_shuffles;
    
    
std::vector<double> sorted_pids = shuffle_pids;

size_t n = sorted_pids.size();
size_t percentile_index = static_cast<size_t>(n * 0.95);

if (percentile_index >= n) percentile_index = n - 1;

std::nth_element(sorted_pids.begin(), 
                sorted_pids.begin() + percentile_index, 
                sorted_pids.end());

double percentile_95 = sorted_pids[percentile_index];
    
    
    if (original_pid < percentile_95) {
        return true;
    } else {
        if (p_value < 0.05) {
            return false; 
        } else {
            return true; 
        }
    }
}


} // namespace dna_utils

#endif 