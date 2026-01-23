// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppParallel)]]
#include "percon2.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <tbb/task_scheduler_init.h>
#include <random>
#include <limits>

using namespace Rcpp;

namespace percon {

std::vector<SequenceData> read_fasta(const std::string& filename) {
  std::vector<SequenceData> sequences;
  std::ifstream file(filename);
  std::string line;
  SequenceData current_seq;
  bool in_sequence = false;
  
  if (!file.is_open()) {
    stop("Cannot open file: " + filename);
  }
  
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    
    if (line[0] == '>') {
      if (in_sequence) {
        sequences.push_back(current_seq);
      }
      
      current_seq.header = line.substr(1);
      current_seq.sequence.clear();
      in_sequence = true;
    } else {
      line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
      std::transform(line.begin(), line.end(), line.begin(), ::toupper);
      current_seq.sequence += line;
    }
  }
  
  if (in_sequence) {
    sequences.push_back(current_seq);
  }
  
  for (size_t i = 0; i < sequences.size(); ++i) {
    sequences[i].index = i;
  }
  
  return sequences;
}

std::unordered_set<std::string> create_kmer_dict(const std::string& sequence, int k) {
  std::unordered_set<std::string> kmer_dict;
  if (sequence.length() < (size_t)k) return kmer_dict;
  
  for (size_t i = 0; i <= sequence.length() - k; ++i) {
    kmer_dict.insert(sequence.substr(i, k));
  }
  return kmer_dict;
}

double calculate_f(const std::unordered_set<std::string>& dict1, 
                   const std::unordered_set<std::string>& dict2) {
  size_t N1 = dict1.size();
  size_t N2 = dict2.size();
  
  if (N1 == 0 || N2 == 0) {
    return std::numeric_limits<double>::infinity();
  }
  
  size_t n_c = 0;
  if (dict1.size() < dict2.size()) {
    for (const auto& kmer : dict1) {
      if (dict2.find(kmer) != dict2.end()) n_c++;
    }
  } else {
    for (const auto& kmer : dict2) {
      if (dict1.find(kmer) != dict1.end()) n_c++;
    }
  }
  
  if (n_c == 0) {
    return std::numeric_limits<double>::infinity();
  }
  
  double value = (N1 + N2) / (2.0 * n_c);
  return std::log(value);
}

std::string generate_random_sequence(int length, const std::string& alphabet) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, alphabet.size() - 1);
  
  std::string random_seq;
  random_seq.reserve(length);
  
  for (int i = 0; i < length; ++i) {
    random_seq += alphabet[dis(gen)];
  }
  return random_seq;
}

double calculate_f_random_single(const std::unordered_set<std::string>& dict, 
                                 int seq_length, int k, int n_iter) {
  double f_sum = 0.0;
  int valid_iterations = 0;
  
  for (int i = 0; i < n_iter; ++i) {
    std::string random_seq = generate_random_sequence(seq_length);
    auto random_dict = create_kmer_dict(random_seq, k);
    double f_val = calculate_f(dict, random_dict);
    
    if (!std::isinf(f_val)) {
      f_sum += f_val;
      valid_iterations++;
    }
  }
  
  if (valid_iterations == 0) {
    return 1.0; 
  }
  
  return f_sum / valid_iterations;
}

PerconWorker::PerconWorker(const std::vector<SequenceData>& sequences,
                           const std::vector<double>& f_random_values,
                           int k,
                           NumericMatrix results,
                           tbb::spin_mutex& mutex)
  : sequences_(sequences), f_random_values_(f_random_values), k_(k), 
    results_(results), mutex_(mutex) {}

void PerconWorker::operator()(std::size_t begin, std::size_t end) {
  for (std::size_t i = begin; i < end; ++i) {
    for (size_t j = i + 1; j < sequences_.size(); ++j) {
      double f_observed = calculate_f(sequences_[i].kmer_dict, sequences_[j].kmer_dict);
      double r_o = std::numeric_limits<double>::infinity();
      
      if (!std::isinf(f_observed)) {
        double avg_f_random = (f_random_values_[i] + f_random_values_[j]) / 2.0;
        r_o = f_observed / avg_f_random;
      }
      
      tbb::spin_mutex::scoped_lock lock(mutex_);
      results_(i, j) = r_o;
      results_(j, i) = r_o; 
    }
    
    results_(i, i) = 0.0;
  }
}

} // namespace percon

// [[Rcpp::export]]
Rcpp::List percon_distance_matrix(const std::string& fasta_file, 
                                  int k = 8, 
                                  int n_iter = 100,
                                  int num_threads = 8,
                                  bool return_mst = 1,
                                bool interspecies = 1) {
  
  if (num_threads > 1) {
    tbb::task_scheduler_init init(num_threads);
  }
  
  Rcout << "Reading FASTA file..." << std::endl;
  std::vector<percon::SequenceData> sequences = percon::read_fasta(fasta_file);
  size_t n_seq = sequences.size();
  Rcout << "Found " << n_seq << " sequences" << std::endl;
  
if (interspecies) {
  Rcout << "Extracting species names from headers..." << std::endl;
  for (auto& seq : sequences) {
    size_t pos = seq.header.find('_');
    if (pos != std::string::npos) {
      seq.species = seq.header.substr(0, pos);
    } else {
      seq.species = seq.header;
    }
  }
}


  Rcout << "Creating k-mer dictionaries..." << std::endl;
  for (auto& seq : sequences) {
    seq.kmer_dict = percon::create_kmer_dict(seq.sequence, k);
  }
  
  Rcout << "Calculating f_random values..." << std::endl;
  std::vector<double> f_random_values(n_seq);
  for (size_t i = 0; i < n_seq; ++i) {
    f_random_values[i] = percon::calculate_f_random_single(sequences[i].kmer_dict, 
                                                           sequences[i].sequence.length(), k, n_iter);
    if (i % 10 == 0) {
      Rcout << "Processed " << i << "/" << n_seq << " sequences" << std::endl;
    }
  }
  
  NumericMatrix results(n_seq, n_seq);
  std::fill(results.begin(), results.end(), NumericVector::get_na());
  
  tbb::spin_mutex mutex;
  
  Rcout << "Calculating distance matrix..." << std::endl;
  percon::PerconWorker worker(sequences, f_random_values, k, results, mutex);
  parallelFor(0, n_seq, worker);
  
StringVector headers(n_seq);
for (size_t i = 0; i < n_seq; ++i) {
  headers[i] = sequences[i].header;
}

rownames(results) = headers;
colnames(results) = headers;

Rcpp::List result = List::create(
  Named("distance_matrix") = results,
  Named("k") = k,
  Named("n_iter") = n_iter,
  Named("sequences") = headers
);

if (interspecies) {
  StringVector species_names(n_seq);
  for (size_t i = 0; i < n_seq; ++i) {
    species_names[i] = sequences[i].species;
  }
  result["species"] = species_names;
}
  
  
  if (return_mst) {
    DataFrame mst_df = build_mst(results);
    
    IntegerVector from_nodes = mst_df["from"];
    IntegerVector to_nodes = mst_df["to"];
    NumericVector weights = mst_df["weight"];
    
    CharacterVector from_names(from_nodes.size());
    CharacterVector to_names(to_nodes.size());
    
    for (int i = 0; i < from_nodes.size(); ++i) {
      from_names[i] = headers[from_nodes[i] - 1];
      to_names[i] = headers[to_nodes[i] - 1];
    }
    
    DataFrame mst_named = DataFrame::create(
      Named("from") = from_names,
      Named("to") = to_names,
      Named("weight") = weights,
      Named("stringsAsFactors") = false
    );
    
    result["mst"] = mst_named;
  }

if (interspecies) {
    Rcout << "Calculating fast PERCON-AMOVA statistics..." << std::endl;
    
    percon::PerconStats stats = percon::calculate_percon_amova_fast(results, sequences);
    
    double p_value = percon::permutation_test_parallel(results, sequences, 100);
    
    result["percon_amova"] = List::create(
        Named("f_statistic") = stats.f_statistic,
        Named("p_value") = p_value,
        Named("total_sum_squares") = stats.total_ss,
        Named("within_group_sum_squares") = stats.within_ss,
        Named("between_group_sum_squares") = stats.between_ss,
        Named("within_group_distance") = stats.within_group_distance,
        Named("between_group_distance") = stats.between_group_distance
    );
    
    Rcout << "PERCON-AMOVA results: F = " << stats.f_statistic 
          << ", p = " << p_value << std::endl;
}

Rcout << "Analyzing graph structure (parallel)..." << std::endl;
percon::Graph filtered_graph = percon::build_filtered_graph(results, 1.0);


std::vector<int> communities = percon::detect_communities_parallel(filtered_graph);
std::vector<double> centrality = percon::calculate_degree_centrality_parallel(filtered_graph);


DataFrame vertex_analysis = percon::create_vertex_analysis(filtered_graph, sequences, communities, centrality);

result["vertex_analysis"] = vertex_analysis;
result["graph_summary"] = List::create(
    Named("num_communities") = *std::max_element(communities.begin(), communities.end()) + 1,
    Named("num_edges") = filtered_graph.edges.size(),
    Named("num_nodes") = filtered_graph.num_nodes
);

Rcout << "Graph analysis complete: " 
      << (*std::max_element(communities.begin(), communities.end()) + 1) 
      << " communities detected" << std::endl;

  return result;
}