// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]
#ifndef PERCON_H
#define PERCON_H

#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_set>
#include <string>
#include <vector>
#include <mutex>
#include <tbb/tbb.h>
#include <random>

using namespace Rcpp;
using namespace RcppParallel;

namespace percon {

struct SequenceData {
  std::string header;
  std::string sequence;
  std::string species;
  size_t index;
  std::unordered_set<std::string> kmer_dict;
};

std::vector<SequenceData> read_fasta(const std::string& filename);

std::unordered_set<std::string> create_kmer_dict(const std::string& sequence, int k);

double calculate_f(const std::unordered_set<std::string>& dict1,
                   const std::unordered_set<std::string>& dict2);

std::string generate_random_sequence(int length, const std::string& alphabet = "ATCG");

double calculate_f_random_single(const std::unordered_set<std::string>& dict,
                                 int seq_length, int k, int n_iter);

class PerconWorker : public Worker {
private:
  const std::vector<SequenceData>& sequences_;
  const std::vector<double>& f_random_values_;
  int k_;
  RMatrix<double> results_;
  tbb::spin_mutex& mutex_;

public:
  PerconWorker(const std::vector<SequenceData>& sequences,
               const std::vector<double>& f_random_values,
               int k,
               NumericMatrix results,
               tbb::spin_mutex& mutex);

  void operator()(std::size_t begin, std::size_t end);
};



struct PerconStats {
    double total_ss;       // Total sum of squares
    double within_ss;      // Within-group sum of squares
    double between_ss;     // Between-group sum of squares
    double f_statistic;    // F-statistic (between/within)
    double within_group_distance;  // Average within-group distance
    double between_group_distance; // Average between-group distance
};

PerconStats calculate_percon_amova(const NumericMatrix& distance_matrix,
                                  const std::vector<SequenceData>& sequences) {
    PerconStats stats = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    std::unordered_map<std::string, std::vector<size_t>> species_groups;
    for (size_t i = 0; i < sequences.size(); ++i) {
        species_groups[sequences[i].species].push_back(i);
    }

    int n_total = sequences.size();
    int n_groups = species_groups.size();

    double total_sum = 0.0;
    int total_pairs = 0;
    for (int i = 0; i < n_total; ++i) {
        for (int j = i + 1; j < n_total; ++j) {
            double dist = distance_matrix(i, j);
            if (!std::isinf(dist)) {
                total_sum += dist * dist;
                total_pairs++;
            }
        }
    }
    stats.total_ss = total_sum;

    double within_sum = 0.0;
    int within_pairs = 0;
    for (const auto& group : species_groups) {
        const auto& indices = group.second;
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = i + 1; j < indices.size(); ++j) {
                double dist = distance_matrix(indices[i], indices[j]);
                if (!std::isinf(dist)) {
                    within_sum += dist * dist;
                    within_pairs++;
                }
            }
        }
    }
    stats.within_ss = within_sum;
    stats.within_group_distance = within_pairs > 0 ? within_sum / within_pairs : 0.0;

    stats.between_ss = stats.total_ss - stats.within_ss;

    double between_sum = 0.0;
    int between_pairs = 0;
    for (const auto& group1 : species_groups) {
        for (const auto& group2 : species_groups) {
            if (group1.first >= group2.first) continue;

            const auto& indices1 = group1.second;
            const auto& indices2 = group2.second;

            for (size_t i : indices1) {
                for (size_t j : indices2) {
                    double dist = distance_matrix(i, j);
                    if (!std::isinf(dist)) {
                        between_sum += dist;
                        between_pairs++;
                    }
                }
            }
        }
    }
    stats.between_group_distance = between_pairs > 0 ? between_sum / between_pairs : 0.0;

    if (stats.within_ss > 0 && n_groups > 1) {
        double ms_between = stats.between_ss / (n_groups - 1);
        double ms_within = stats.within_ss / (n_total - n_groups);
        stats.f_statistic = ms_between / ms_within;
    }

    return stats;
}

double permutation_test(const NumericMatrix& distance_matrix,
                       const std::vector<SequenceData>& sequences,
                       int n_permutations = 1000) {
    PerconStats observed_stats = calculate_percon_amova(distance_matrix, sequences);
    double observed_f = observed_stats.f_statistic;

    if (observed_f <= 0) return 1.0;

    int count_extreme = 0;
    std::vector<SequenceData> permuted_sequences = sequences;
    std::vector<std::string> original_species;

    for (const auto& seq : sequences) {
        original_species.push_back(seq.species);
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int perm = 0; perm < n_permutations; ++perm) {
        std::vector<std::string> shuffled_species = original_species;
        std::shuffle(shuffled_species.begin(), shuffled_species.end(), gen);

        for (size_t i = 0; i < permuted_sequences.size(); ++i) {
            permuted_sequences[i].species = shuffled_species[i];
        }

        PerconStats perm_stats = calculate_percon_amova(distance_matrix, permuted_sequences);

        if (perm_stats.f_statistic >= observed_f) {
            count_extreme++;
        }
    }

    return (double)(count_extreme + 1) / (n_permutations + 1);
}

struct GraphEdge {
    int from;
    int to;
    double weight;

    GraphEdge(int f, int t, double w) : from(f), to(t), weight(w) {}
};

struct Graph {
    std::vector<GraphEdge> edges;
    int num_nodes;
};

Graph build_filtered_graph(const NumericMatrix& distance_matrix, double threshold = 1.0) {
    Graph graph;
    graph.num_nodes = distance_matrix.nrow();

    for (int i = 0; i < graph.num_nodes; ++i) {
        for (int j = i + 1; j < graph.num_nodes; ++j) {
            double ro = distance_matrix(i, j);

            if (!std::isinf(ro) && ro < threshold) {
                graph.edges.emplace_back(i, j, ro);
            }
        }
    }

    return graph;
}

std::vector<int> detect_communities(const Graph& graph, int max_iterations = 100) {
    std::vector<int> communities(graph.num_nodes);

    for (int i = 0; i < graph.num_nodes; ++i) {
        communities[i] = i;
    }

    bool changed = true;
    int iterations = 0;

    while (changed && iterations < max_iterations) {
        changed = false;

        std::vector<int> node_order(graph.num_nodes);
        for (int i = 0; i < graph.num_nodes; ++i) node_order[i] = i;
        std::random_shuffle(node_order.begin(), node_order.end());

        for (int node : node_order) {
            std::unordered_map<int, double> community_weights;

            for (const auto& edge : graph.edges) {
                int neighbor = -1;
                if (edge.from == node) neighbor = edge.to;
                if (edge.to == node) neighbor = edge.from;

                if (neighbor != -1) {
                    int neighbor_community = communities[neighbor];
                    community_weights[neighbor_community] += 1.0 / edge.weight; // Чем меньше вес, тем сильнее связь
                }
            }

            if (!community_weights.empty()) {
                int best_community = std::max_element(
                    community_weights.begin(), community_weights.end(),
                    [](const auto& a, const auto& b) { return a.second < b.second; }
                )->first;

                if (best_community != communities[node]) {
                    communities[node] = best_community;
                    changed = true;
                }
            }
        }
        iterations++;
    }

    return communities;
}

std::vector<double> calculate_degree_centrality(const Graph& graph) {
    std::vector<double> centrality(graph.num_nodes, 0.0);
    std::vector<int> degree(graph.num_nodes, 0);

    for (const auto& edge : graph.edges) {
        degree[edge.from]++;
        degree[edge.to]++;
    }

    int max_degree = *std::max_element(degree.begin(), degree.end());
    if (max_degree > 0) {
        for (int i = 0; i < graph.num_nodes; ++i) {
            centrality[i] = static_cast<double>(degree[i]) / max_degree;
        }
    }

    return centrality;
}

DataFrame create_vertex_analysis(const Graph& graph,
                                const std::vector<SequenceData>& sequences,
                                const std::vector<int>& communities,
                                const std::vector<double>& centrality) {

    StringVector vertex_names(graph.num_nodes);
    StringVector species_names(graph.num_nodes);
    IntegerVector community_ids(graph.num_nodes);
    NumericVector centrality_scores(graph.num_nodes);

    for (int i = 0; i < graph.num_nodes; ++i) {
        vertex_names[i] = sequences[i].header;
        species_names[i] = sequences[i].species;
        community_ids[i] = communities[i] + 1; // 1-based
        centrality_scores[i] = centrality[i];
    }

    return DataFrame::create(
        Named("Vertex") = vertex_names,
        Named("Species") = species_names,
        Named("Community") = community_ids,
        Named("Centrality") = centrality_scores,
        Named("stringsAsFactors") = false
    );
}


PerconStats calculate_percon_amova_fast(const NumericMatrix& distance_matrix,
                                       const std::vector<SequenceData>& sequences) {
    PerconStats stats = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    std::unordered_map<std::string, std::vector<size_t>> species_groups;
    for (size_t i = 0; i < sequences.size(); ++i) {
        species_groups[sequences[i].species].push_back(i);
    }

    int n_total = sequences.size();
    int n_groups = species_groups.size();

    double total_sum = 0.0;
    int total_count = 0;
    for (int i = 0; i < n_total; ++i) {
        for (int j = i + 1; j < n_total; ++j) {
            double dist = distance_matrix(i, j);
            if (!std::isinf(dist)) {
                total_sum += dist * dist;
                total_count++;
            }
        }
    }
    stats.total_ss = total_sum;

    double within_sum = 0.0;
    int within_count = 0;
    for (const auto& group : species_groups) {
        const auto& indices = group.second;
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = i + 1; j < indices.size(); ++j) {
                double dist = distance_matrix(indices[i], indices[j]);
                if (!std::isinf(dist)) {
                    within_sum += dist * dist;
                    within_count++;
                }
            }
        }
    }
    stats.within_ss = within_sum;
    stats.between_ss = total_sum - within_sum;

    stats.within_group_distance = within_count > 0 ? within_sum / within_count : 0.0;
    stats.between_group_distance = (total_count - within_count) > 0 ?
                                  stats.between_ss / (total_count - within_count) : 0.0;

    if (stats.within_ss > 0 && n_groups > 1) {
        double ms_between = stats.between_ss / (n_groups - 1);
        double ms_within = stats.within_ss / (n_total - n_groups);
        stats.f_statistic = ms_between / ms_within;
    }

    return stats;
}


double permutation_test_parallel(const NumericMatrix& distance_matrix,
                                const std::vector<SequenceData>& sequences,
                                int n_permutations = 100) {
    PerconStats observed_stats = calculate_percon_amova_fast(distance_matrix, sequences);
    double observed_f = observed_stats.f_statistic;

    if (observed_f <= 0) return 1.0;

    std::vector<std::string> original_species;
    for (const auto& seq : sequences) {
        original_species.push_back(seq.species);
    }

    tbb::atomic<int> count_extreme;
    count_extreme = 0;

    class PermutationWorker {
    public:
        PermutationWorker(const NumericMatrix& dist_matrix,
                         const std::vector<SequenceData>& seqs,
                         const std::vector<std::string>& orig_species,
                         double observed_f_val,
                         tbb::atomic<int>& counter)
            : distance_matrix_(dist_matrix), sequences_(seqs),
              original_species_(orig_species), observed_f_(observed_f_val),
              count_extreme_(counter) {}

        void operator()(const tbb::blocked_range<int>& range) const {
            static thread_local std::mt19937 gen(std::random_device{}());

            for (int perm = range.begin(); perm != range.end(); ++perm) {
                std::vector<std::string> shuffled_species = original_species_;
                std::shuffle(shuffled_species.begin(), shuffled_species.end(), gen);

                std::vector<SequenceData> temp_sequences = sequences_;
                for (size_t i = 0; i < temp_sequences.size(); ++i) {
                    temp_sequences[i].species = shuffled_species[i];
                }

                PerconStats perm_stats = calculate_percon_amova_fast(distance_matrix_, temp_sequences);

                if (perm_stats.f_statistic >= observed_f_) {
                    count_extreme_++;
                }
            }
        }

    private:
        const NumericMatrix& distance_matrix_;
        const std::vector<SequenceData>& sequences_;
        const std::vector<std::string>& original_species_;
        double observed_f_;
        tbb::atomic<int>& count_extreme_;
    };

    PermutationWorker worker(distance_matrix, sequences, original_species, observed_f, count_extreme);
    tbb::parallel_for(tbb::blocked_range<int>(0, n_permutations), worker);

    return (double)(count_extreme + 1) / (n_permutations + 1);
}






std::vector<int> detect_communities_parallel(const Graph& graph, double min_strength = 0.6, int max_iterations = 50) {
    std::vector<int> communities(graph.num_nodes);

    for (int i = 0; i < graph.num_nodes; ++i) {
        communities[i] = i;  // Узел i → сообщество i
    }

    std::vector<std::vector<int>> neighbors(graph.num_nodes);
    std::vector<std::vector<double>> neighbor_weights(graph.num_nodes);

    for (const auto& edge : graph.edges) {
        double strength = 1.0 / (edge.weight + 0.001);
       if (strength >= min_strength) {
            neighbors[edge.from].push_back(edge.to);
            neighbor_weights[edge.from].push_back(strength);

            neighbors[edge.to].push_back(edge.from);
            neighbor_weights[edge.to].push_back(strength);
        }
    }

    bool changed = true;
    int iterations = 0;

    while (changed && iterations < max_iterations) {
        changed = false;
        tbb::atomic<bool> global_changed;
        global_changed = false;

        std::vector<int> node_order(graph.num_nodes);
        for (int i = 0; i < graph.num_nodes; ++i) node_order[i] = i;
        std::random_shuffle(node_order.begin(), node_order.end());

        class CommunityWorker {
        public:
            CommunityWorker(const std::vector<std::vector<int>>& neigh,
                          const std::vector<std::vector<double>>& weights,
                          const std::vector<int>& order,
                          std::vector<int>& comm,
                          tbb::atomic<bool>& changed_flag)
                : neighbors_(neigh), neighbor_weights_(weights),
                  node_order_(order), communities_(comm), global_changed_(changed_flag) {}

            void operator()(const tbb::blocked_range<int>& range) const {
                for (int idx = range.begin(); idx != range.end(); ++idx) {
                    int node = node_order_[idx];

                    if (neighbors_[node].empty()) continue;

                    std::unordered_map<int, double> community_weights;

                    for (size_t i = 0; i < neighbors_[node].size(); ++i) {
                        int neighbor = neighbors_[node][i];
                        double strength = neighbor_weights_[node][i];
                        int neighbor_community = communities_[neighbor];
                        community_weights[neighbor_community] += strength;
                    }

                    if (!community_weights.empty()) {
                        auto best_community = std::max_element(
                            community_weights.begin(), community_weights.end(),
                            [](const auto& a, const auto& b) { return a.second < b.second; }
                        );

                        if (best_community->first != communities_[node]) {
                            communities_[node] = best_community->first;
                            global_changed_ = true;
                        }
                    }
                }
            }

        private:
            const std::vector<std::vector<int>>& neighbors_;
            const std::vector<std::vector<double>>& neighbor_weights_;
            const std::vector<int>& node_order_;
            std::vector<int>& communities_;
            tbb::atomic<bool>& global_changed_;
        };

        CommunityWorker worker(neighbors, neighbor_weights, node_order, communities, global_changed);
        tbb::parallel_for(tbb::blocked_range<int>(0, graph.num_nodes), worker);

        changed = global_changed;
        iterations++;

        Rcout << "Iteration " << iterations << ": "
              << (changed ? "changes occurred" : "converged") << std::endl;

        if (!changed) break;
    }

    std::unordered_map<int, int> community_remap;
    int new_id = 0;
    for (int& comm : communities) {
        if (community_remap.find(comm) == community_remap.end()) {
            community_remap[comm] = new_id++;
        }
        comm = community_remap[comm];
    }

    std::vector<int> community_sizes(new_id, 0);
    for (int comm : communities) {
        community_sizes[comm]++;
    }

    Rcout << "Detected " << new_id << " communities" << std::endl;
    Rcout << "Community size distribution: ";
    for (int size : community_sizes) {
        Rcout << size << " ";
    }
    Rcout << std::endl;

    return communities;
}




std::vector<double> calculate_degree_centrality_parallel(const Graph& graph) {
    std::vector<double> centrality(graph.num_nodes, 0.0);
    std::vector<int> degree(graph.num_nodes, 0);

    class DegreeWorker {
    public:
        DegreeWorker(const Graph& g, std::vector<int>& deg) : graph_(g), degree_(deg) {}

        void operator()(const tbb::blocked_range<size_t>& range) const {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                const auto& edge = graph_.edges[i];
                degree_[edge.from]++;
                degree_[edge.to]++;
            }
        }

    private:
        const Graph& graph_;
        std::vector<int>& degree_;
    };

    DegreeWorker degree_worker(graph, degree);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, graph.edges.size()), degree_worker);

    int max_degree = 0;
    for (int d : degree) {
        if (d > max_degree) max_degree = d;
    }

    class CentralityWorker {
    public:
        CentralityWorker(std::vector<double>& cent, const std::vector<int>& deg, int max_deg)
            : centrality_(cent), degree_(deg), max_degree_(max_deg) {}

        void operator()(const tbb::blocked_range<int>& range) const {
            for (int i = range.begin(); i != range.end(); ++i) {
                centrality_[i] = static_cast<double>(degree_[i]) / max_degree_;
            }
        }

    private:
        std::vector<double>& centrality_;
        const std::vector<int>& degree_;
        int max_degree_;
    };

    if (max_degree > 0) {
        CentralityWorker centrality_worker(centrality, degree, max_degree);
        tbb::parallel_for(tbb::blocked_range<int>(0, graph.num_nodes), centrality_worker);
    }

    return centrality;
}




} // namespace percon

Rcpp::List percon_distance_matrix(const std::string& fasta_file,
                                  int k = 8,
                                  int n_iter = 100,
                                  int num_threads = 1,
                                  bool return_mst = true);

struct MSTEdge {
  int from;
  int to;
  double weight;

  bool operator>(const MSTEdge& other) const {
    return weight > other.weight;
  }
};


Rcpp::DataFrame build_mst(const Rcpp::NumericMatrix& dist_matrix) {
  int n = dist_matrix.rows();
  std::vector<bool> in_mst(n, false);
  std::vector<MSTEdge> mst_edges;
  std::priority_queue<MSTEdge, std::vector<MSTEdge>, std::greater<MSTEdge>> pq;

  in_mst[0] = true;
  for (int j = 1; j < n; ++j) {
    if (!std::isinf(dist_matrix(0, j))) {
      pq.push({0, j, dist_matrix(0, j)});
    }
  }

  while (!pq.empty() && mst_edges.size() < n - 1) {
    MSTEdge current = pq.top();
    pq.pop();

    if (in_mst[current.to]) continue;

    in_mst[current.to] = true;
    mst_edges.push_back(current);

    for (int j = 0; j < n; ++j) {
      if (!in_mst[j] && !std::isinf(dist_matrix(current.to, j))) {
        pq.push({current.to, j, dist_matrix(current.to, j)});
      }
    }
  }

  Rcpp::IntegerVector from(mst_edges.size()), to(mst_edges.size());
  Rcpp::NumericVector weight(mst_edges.size());

  for (size_t i = 0; i < mst_edges.size(); ++i) {
    from[i] = mst_edges[i].from + 1; //  1-based
    to[i] = mst_edges[i].to + 1;
    weight[i] = mst_edges[i].weight;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("from") = from,
    Rcpp::Named("to") = to,
    Rcpp::Named("weight") = weight
  );
}

#endif
