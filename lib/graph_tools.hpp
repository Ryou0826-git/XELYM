// graph_tools.hpp

#ifndef GRAPH_TOOLS_HPP
#define GRAPH_TOOLS_HPP

#include <igraph.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>

#include <psf_tools.hpp>

using namespace std;

struct s_graph {
    igraph_t        graph;
    vector<int>     node_names;
    igraph_vector_t weights;
    igraph_matrix_t spl;
};

struct s_gn {
     igraph_t            graph;
     igraph_vector_int_t membership;
     igraph_integer_t    ncomp;
     igraph_real_t       modularity;
     vector<int>         node_names;
};

bool edge_exists(const igraph_t& graph, int from, int to);

int find_node_by_label(const s_graph& G, int label_to_find);

igraph_t load_graph_from_gml(const std::string& filename);
void export_graph_to_gml(const igraph_t* graph,
                         const std::string& filename);

void print_graph_summary(const igraph_t& graph);

void deep_copy_graph(s_graph& G, const s_graph& G_bond);

void DestroyGN(s_gn& gn);
void DestroyGRAPH(s_graph& G);

///////////////////////////////////////////////////////////////////////////////

igraph_matrix_t CalcSPL(igraph_t& graph, igraph_vector_t& weights);
vector<vector<double>>
CalcSPLaveCommunity(s_graph& G,
                    const igraph_vector_int_t& membership);

pair<int, double> get_max_betweenness_edge(const igraph_t& graph);

s_gn girvan_newman(const s_graph& G,
                   const int& nedge);

s_graph GenerateBond(const s_psf& psf,
                     const vector<int>& labellist);

s_graph GenerateCommunity(const s_graph& G_bond,
                          const vector<vector<double>>& coord,
                          const double& cutoff);

s_graph GenerateCommunityEdge(const vector<vector<int>>& resinfo,
                              const vector<int>& residinfo,
                              const vector<vector<double>>& coord,
                              const double& cutoff,
                              const vector<vector<double>>& wij);

s_graph GenerateCommunityLen(const vector<vector<int>>& resinfo,
                             const vector<int>& residinfo,
                             const vector<vector<double>>& coord,
                             const double& cutoff);

//
// for NMI, ARI
vector<vector<int>> confusion_matrix(const vector<int>& memb_refs,
                                     const vector<int>& memb);

double normalized_mutual_information(const vector<int>& memb_refs,
                                     const vector<int>& memb);

double adjusted_rand_index(const vector<int>& memb_refs,
                           const vector<int>& memb);

pair<int, double> ReturnCommNMI(const vector<vector<int>>& refs_memb,
                                const vector<int>& memb);

#endif
