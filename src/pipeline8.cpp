// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]


#include <Rcpp.h>
#include <RcppParallel.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "dna_utils.h"
#include "CoresScan.h"
#include "Extend.h"
#include "repsplit.h"
#include "compfilt.h"
#include "refinement.h"
#include "smartfilter.h"
#include "spacerFilter.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace dna_utils;

std::string format_time(int seconds) {
  int hours = seconds / 3600;
  int minutes = (seconds % 3600) / 60;
  int secs = seconds % 60;
  
  if (hours > 0) {
    return std::to_string(hours) + "h " + std::to_string(minutes) + "m " + std::to_string(secs) + "s";
  } else if (minutes > 0) {
    return std::to_string(minutes) + "m " + std::to_string(secs) + "s";
  } else {
    return std::to_string(secs) + "s";
  }
}

struct PipelineResult {
  std::string consensus;
  int dr_length;
  std::vector<int> positions;
  std::vector<std::string> sequences;
  std::string seqname; 
  
  PipelineResult() : dr_length(0) {};
  
  PipelineResult(const CoreResult& core_result, const std::string& seqname)
    : consensus(core_result.consensus),
      dr_length(core_result.dr_length),
      positions(core_result.positions),
      sequences(core_result.sequences),
      seqname(seqname) {}
};

inline std::pair<std::string, std::string> read_first_sequence_info(const std::string& fasta_path) {
  std::ifstream file(fasta_path);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open FASTA file: " + fasta_path);
  }
  
  std::string line;
  std::string seqname;
  std::string first_chunk;
  bool in_sequence = false;
  int chunk_size = 0;
  const int max_chunk_size = 10000; 
  
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    
    if (line[0] == '>') {
      if (!seqname.empty()) break; 
      
      std::istringstream header_line(line.substr(1));
      header_line >> seqname;
    } else if (!seqname.empty()) {
      first_chunk += line;
      chunk_size += line.size();
      if (chunk_size >= max_chunk_size) break;
    }
  }
  
  file.close();
  
  if (seqname.empty()) {
    throw std::runtime_error("No sequences found in FASTA file");
  }
  
  return {seqname, first_chunk};
}

std::vector<std::pair<std::string, std::string>> read_fasta_by_chunks(
    const std::string& fasta_path, 
    int chunk_size,
    int overlap) {
  
  std::ifstream file(fasta_path);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open FASTA file: " + fasta_path);
  }
  
  std::vector<std::pair<std::string, std::string>> chunks;
  std::string line;
  std::string current_seqname;
  std::string current_chunk;
  int current_size = 0;
  
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    
    if (line[0] == '>') {
      if (!current_chunk.empty()) {
        chunks.emplace_back(current_seqname, current_chunk);
        current_chunk.clear();
        current_size = 0;
      }
      
      std::istringstream header_line(line.substr(1));
      header_line >> current_seqname;
    } else {
      if (current_size >= chunk_size) {
        std::string chunk_with_overlap = current_chunk.substr(
          current_chunk.size() - overlap); 
        
        chunks.emplace_back(current_seqname, current_chunk);
        current_chunk = chunk_with_overlap; 
        current_size = chunk_with_overlap.size();
      }
      
      current_chunk += line;
      current_size += line.size();
    }
  }
  
  if (!current_chunk.empty()) {
    chunks.emplace_back(current_seqname, current_chunk);
  }
  
  file.close();
  return chunks;
}

inline void adjust_positions(std::vector<PipelineResult>& results, 
                             const std::vector<int>& chunk_offsets,
                             const std::vector<std::string>& chunk_seqnames) {
  for (size_t i = 0; i < results.size(); ++i) {
    for (auto& pos : results[i].positions) {
      for (size_t j = 0; j < chunk_offsets.size(); ++j) {
        if (pos >= chunk_offsets[j] && 
            (j == chunk_offsets.size() - 1 || pos < chunk_offsets[j + 1])) {
          pos -= chunk_offsets[j];
          break;
        }
      }
    }
    
    for (size_t j = 0; j < chunk_offsets.size(); ++j) {
      if (results[i].positions[0] >= chunk_offsets[j] && 
          (j == chunk_offsets.size() - 1 || results[i].positions[0] < chunk_offsets[j + 1])) {
        results[i].seqname = chunk_seqnames[j];
        break;
      }
    }
  }
}




List convert_to_r_list(const std::vector<PipelineResult>& results) {
  List r_results(results.size());
  
  for (size_t i = 0; i < results.size(); ++i) {
    List core_list;
    core_list["consensus"] = results[i].consensus;
    core_list["dr_length"] = results[i].dr_length;
    
    // Конвертируем позиции в 1-based для R!
    std::vector<int> positions_1based = results[i].positions;
    for (auto& pos : positions_1based) {
        pos += 1; // 0-based → 1-based
    }
    core_list["positions"] = wrap(positions_1based);
    
    core_list["sequences"] = wrap(results[i].sequences);
    core_list["seqname"] = results[i].seqname;
    r_results[i] = core_list;
  }
  
  return r_results;
}



void log_final_results(const std::vector<PipelineResult>& results, 
                      const std::string& full_genome,
                      const std::string& filename) {
    std::ofstream logfile(filename);
    
    logfile << "Element\tStart\tStop\tCassette_id\tSequence\tSeqname\n";
    
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& result = results[i];
        int cassette_id = i + 1;
        
        if (result.positions.size() < 2) {
            continue;
        }
        
        std::vector<int> sorted_positions = result.positions;
        std::sort(sorted_positions.begin(), sorted_positions.end());
        
        auto spacers = extract_spacers(full_genome, sorted_positions, result.dr_length);
        
        int cassette_start = sorted_positions.front();
        int cassette_end = sorted_positions.back() + result.dr_length;
        
        logfile << "Cassette\t" << cassette_start << "\t" << cassette_end 
                << "\t" << cassette_id << "\t" << "FULL_CASSETTE" << "\t" << result.seqname << "\n";
        
        for (size_t j = 0; j < sorted_positions.size(); ++j) {
            int dr_start = sorted_positions[j];
            int dr_end = dr_start + result.dr_length;
            
            logfile << "DR\t" << dr_start << "\t" << dr_end 
                    << "\t" << cassette_id << "\t" << result.consensus << "\t" << result.seqname << "\n";
            
            if (j < spacers.size()) {
                int sp_start = dr_end;
                int sp_end = sorted_positions[j + 1];
                
                logfile << "SP\t" << sp_start << "\t" << sp_end 
                        << "\t" << cassette_id << "\t" << spacers[j] << "\t" << result.seqname << "\n";
            }
        }
        
        logfile << "---\n";
    }
    
    logfile.close();
}

std::string read_complete_genome(const std::string& fasta_path) {
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open genome FASTA: " + fasta_path);
    }
    
    std::string genome_sequence;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue; 
        
        genome_sequence += line;
    }
    
    file.close();
    Rcout << "Complete genome loaded: " << genome_sequence.size() << " bp" << std::endl;
    return genome_sequence;
}


void check_stage_leaks(const std::vector<CoreResult>& results, const std::string& stage_name) {
    for (const auto& result : results) {
        if (result.dr_length <= 0 || result.dr_length > 100) {
            std::ofstream leak_trap("PIPELINE_LEAKS.txt", std::ios_base::app);
            leak_trap << "LEAK at stage: " << stage_name << std::endl;
            leak_trap << "dr_length: " << result.dr_length << std::endl;
            leak_trap << "consensus: " << result.consensus << std::endl;
            leak_trap << "consensus_length: " << result.consensus.length() << std::endl;
            leak_trap << "positions_count: " << result.positions.size() << std::endl;
            leak_trap << "---" << std::endl;
        }
    }
}



// [[Rcpp::export]]
List dna_pipeline_master(
    const std::string& fasta_path,
    int chunk_size = 1000000,
    int overlap = 1000,
    int k_core = 10,
    int min_repeats = 3,
    int window_size = 2000,
    int step_size = 1000,
    int num_threads = 4,
    int min_matches = 8,
    int max_search_offset = 1500,
    int max_extension_length = 55,
    int entropy_window_size = 10,
    int max_k = 8,
    double dr_threshold = 0.5,
    double spacer_threshold = 0.6,
    double min_similarity = 0.8,
    double similarity_threshold = 0.9,
    double spacer_similarity_threshold = 0.3,
    bool heuristics = true) {
  
  auto chunks = read_fasta_by_chunks(fasta_path, chunk_size, overlap);
  int total_chunks = chunks.size();
  
  std::vector<PipelineResult> all_results;
  std::vector<int> chunk_offsets;
  std::vector<std::string> chunk_seqnames;
  
  int current_offset = 0;
  std::string current_sequence_name = "";
  
  auto start_time = std::chrono::steady_clock::now();
  


  for (int chunk_idx = 0; chunk_idx < total_chunks; ++chunk_idx) {
    const auto& chunk = chunks[chunk_idx];
    const auto& seqname = chunk.first;
    const auto& sequence = chunk.second;
    
    if (seqname != current_sequence_name) {
        current_offset = 0;
        current_sequence_name = seqname;
        Rcout << "NEW SEQUENCE: " << seqname << " - reset offset to 0" << std::endl;
    }
    
    auto current_time = std::chrono::steady_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
    double speed = chunk_idx > 0 ? elapsed_seconds / static_cast<double>(chunk_idx) : 0;
    int estimated_seconds_left = static_cast<int>(speed * (total_chunks - chunk_idx));
    
    Rcout << "Working on chunk " << (chunk_idx + 1) << " of " << total_chunks 
          << " | Elapsed: " << format_time(elapsed_seconds)
          << " | Estimated left: ~" << format_time(estimated_seconds_left)
          << " | Sequence: " << seqname << std::endl;
    
    chunk_offsets.push_back(current_offset);
    chunk_seqnames.push_back(seqname);
    
    auto chunk_results = deep_cores_search_internal(
      sequence, k_core, min_repeats, window_size, step_size, 
      num_threads, min_matches, max_search_offset, max_extension_length, heuristics);

    
    auto extended_results = extend_cores_internal(chunk_results, sequence, max_extension_length,
                                                  entropy_window_size, num_threads);


    auto split_results = split_cassettes_internal(extended_results);



    auto filtered_results = complexity_filter_internal(split_results, sequence, max_k, dr_threshold, spacer_threshold, num_threads);


    auto refined_results = refine_cassettes_internal(filtered_results, sequence, min_similarity, spacer_similarity_threshold, num_threads);


    auto sfinal_results = smart_filter_internal(refined_results, sequence, similarity_threshold, num_threads);


    auto final_results = spacer_filter_internal(sfinal_results, sequence);


    for (auto& core_result : final_results) {
      PipelineResult pipeline_result(core_result, seqname);
      
      for (auto& pos : pipeline_result.positions) {
        pos += current_offset;
      }
      
      all_results.push_back(std::move(pipeline_result));
    }
    
if (chunk_idx < total_chunks - 1) {
  current_offset += (sequence.size() - overlap); 
} else {
  current_offset += sequence.size();
}
  }
  
 
  return convert_to_r_list(all_results);
}