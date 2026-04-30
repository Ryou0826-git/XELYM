// graph_tools.cpp

#include <iomanip> 
#include <graph_tools.hpp>

//---------------------------------------------------------------------------//
bool edge_exists(const igraph_t& graph, int from, int to) {
    igraph_integer_t eid;
    igraph_get_eid(&graph, &eid, from, to, IGRAPH_UNDIRECTED, false);
    //
    if (eid != -1) {
        return true;
    } else {
        return false;
    }
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
int find_node_by_label(const s_graph& G, int label_to_find) {
    int nnode = G.node_names.size();
    //
    for (int i = 0; i < nnode; i++) {
        int name = G.node_names[i];
        if (name == label_to_find) {
            return i;
        }
    }

    return -1;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
igraph_t load_graph_from_gml(const string& filename) {
    igraph_t graph;
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);

    FILE* file = fopen(filename.c_str(), "r");
    if (!file) {
        cerr << "Failed to open file " << filename << endl;
        return graph;
    }

    if (igraph_read_graph_gml(&graph, file) != IGRAPH_SUCCESS) {
        cerr << "Failed to read graph from GML" << endl;
        fclose(file);
        return graph;
    }

    fclose(file);
    return graph;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void export_graph_to_gml(const igraph_t* graph, const string& filename) {
    FILE* file = fopen(filename.c_str(), "w");
    if (!file) {
        cerr << "Failed to open file " << filename << endl;
        return;
    }

    if (igraph_write_graph_gml(graph, file, IGRAPH_WRITE_GML_DEFAULT_SW, 
                               nullptr, nullptr) != IGRAPH_SUCCESS) {
        cerr << "Failed to write graph to GML" << endl;
    } else {
        cout << "Graph written to " << filename << endl;
    }

    fclose(file);
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void print_graph_summary(const igraph_t& graph) {
    //
    int num_vertices = igraph_vcount(&graph);
    int num_edges    = igraph_ecount(&graph);

    double density   = (2.0 * num_edges) / (num_vertices * (num_vertices - 1));

    igraph_vector_int_t degree;
    igraph_vector_int_init(&degree, 0);
    igraph_degree(&graph, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

    int max_degree = VECTOR(degree)[0];
    int min_degree = VECTOR(degree)[0];
    for (int i = 1; i < num_vertices; ++i) {
        if (VECTOR(degree)[i] > max_degree) {
            max_degree = VECTOR(degree)[i];
        }
        if (VECTOR(degree)[i] < min_degree) {
            min_degree = VECTOR(degree)[i];
        }
    }

    igraph_real_t diameter;
    igraph_diameter(&graph, &diameter, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, true);

    cout << endl;
    cout << "Number of vertices: " << num_vertices << endl;
    cout << "Number of edges: "    << num_edges    << endl;
    cout << "Density: "            << density      << endl;
    cout << "Max degree: "         << max_degree   << endl;
    cout << "Min degree: "         << min_degree   << endl;
    cout << "Diameter: "           << diameter     << endl;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void DestroyGN(s_gn& gn) {
    igraph_vector_int_destroy(&gn.membership);
    igraph_destroy(&gn.graph);
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void DestroyGRAPH(s_graph& G) {
    igraph_vector_destroy(&G.weights);
    igraph_matrix_destroy(&G.spl);
    igraph_destroy(&G.graph);
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void deep_copy_graph(s_graph& G, const s_graph& G_bond) {
    //
    igraph_copy(&G.graph, &G_bond.graph);
    G.node_names = G_bond.node_names;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
igraph_matrix_t CalcSPL(igraph_t& graph,
                        igraph_vector_t& weights) {
    igraph_matrix_t res;
    //
    igraph_matrix_init(&res, igraph_vcount(&graph), igraph_vcount(&graph));
    igraph_shortest_paths_dijkstra(&graph, &res, igraph_vss_all(), igraph_vss_all(),
                                   &weights, IGRAPH_ALL);

    return res;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
vector<vector<double>>
CalcSPLaveCommunity(s_graph& G,
                    const igraph_vector_int_t& membership) {
    //
    int ncomp = -1;
    int lmemb = igraph_vector_int_size(&membership);
    for (int i = 0; i < lmemb; ++i) {
        if (ncomp < VECTOR(membership)[i]) {
            ncomp = VECTOR(membership)[i];
        }
    }
    ncomp += 1;
    //

    G.spl = CalcSPL(G.graph, G.weights);
    vector<vector<double>> splave(ncomp, vector<double>(ncomp, 0.0));
    vector<vector<double>> count(ncomp, vector<double>(ncomp, 0.0));

    for (int i = 0; i < lmemb; ++i) {
        int member_i = VECTOR(membership)[i];
        for (int j = 0; j < lmemb; ++j) {
            int member_j = VECTOR(membership)[j];
            double spl = MATRIX(G.spl, i, j);
            splave[member_i][member_j] += spl;
            count[member_i][member_j]  += 1;
        }
    }

    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            splave[i][j] /= count[i][j];
        }
    }

    return splave;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
pair<int, double> 
get_max_betweenness_edge(const igraph_t& graph,
                         const igraph_vector_t& weights) {
    igraph_vector_t betweenness;
    igraph_vector_init(&betweenness, 0);

    //cout << VECTOR(weights)[2] << endl;

    igraph_edge_betweenness(&graph, &betweenness, IGRAPH_UNDIRECTED, &weights);

    igraph_real_t max_value = -1;
    int max_index = -1;
    for (int i = 0; i < igraph_ecount(&graph); ++i) {
        if (VECTOR(betweenness)[i] > max_value) {
            max_value = VECTOR(betweenness)[i];
            max_index = i;
        }
    }
    
    //cout << fixed << setprecision(10);
    //for (int i = 0; i < igraph_ecount(&graph); i++) {
    //    cout << "Edge " << i << ": " << VECTOR(betweenness)[i] << endl;
    //}

    igraph_vector_destroy(&betweenness);

    return {max_index, max_value};
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_gn girvan_newman(const s_graph& G,
                   const int& ncomp) {
    //
    s_gn gn;
    //
    igraph_vector_int_t tmp_membership;
    igraph_vector_int_init(&tmp_membership, 0);

    igraph_integer_t num_components = -1;
    igraph_real_t modularity        = -1;
    igraph_real_t best_modularity   = -1;

    igraph_t        tmp_graph;
    igraph_vector_t tmp_weights;

    igraph_copy(&tmp_graph, &G.graph);
    igraph_vector_copy(&tmp_weights, &G.weights);

    gn.node_names = G.node_names;
 
    cout << endl;
    while (num_components < ncomp) {
        //
        auto [edge_index, edge_value] = get_max_betweenness_edge(tmp_graph, tmp_weights);
        if (edge_index == -1) break;

        igraph_es_t es;
        igraph_es_1(&es, edge_index);
        igraph_delete_edges(&tmp_graph, es);
        igraph_es_destroy(&es);
        igraph_vector_remove(&tmp_weights, edge_index);

        cout << "deleted " << edge_index << " edge" << endl;

        igraph_connected_components(&tmp_graph, &tmp_membership, NULL,
                                    &num_components, IGRAPH_WEAK);

        igraph_modularity(&tmp_graph, &tmp_membership, &tmp_weights,
                          1.0, IGRAPH_UNDIRECTED, &modularity);
        //
        if (modularity > best_modularity) {
            best_modularity = modularity;
        }
        //
        cout << "Current modularity: " << modularity
             << ", Components: "       << num_components << endl;
        cout << "edge betweenness: "   << edge_value     << endl;
    }

    //
    // Copy graph info.
    gn.modularity = best_modularity;
    gn.ncomp      = num_components;
    igraph_copy(&gn.graph, &tmp_graph);
    igraph_vector_int_copy(&gn.membership, &tmp_membership);

    cout << "Best modularity: " << best_modularity << endl;

    vector<vector<int>> ClusterSel(ncomp);
    for (int i = 0; i < igraph_vector_int_size(&gn.membership); ++i) {
        int sel = VECTOR(gn.membership)[i];
        ClusterSel[sel].push_back(gn.node_names[i]);
    }

    for (int i = 0; i < ncomp; ++i) {
        int cnum = ClusterSel[i].size();
        cout << "Cluster " << i+1 << " : ";
        for (int j = 0; j < cnum; ++j) {
            cout << ClusterSel[i][j] << " ";
        }
        cout << endl;
    }

    return gn;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_graph GenerateBond(const s_psf& psf,
                     const vector<int>& labellist) {
    s_graph G;
    //
    int nnode = labellist.size();
    igraph_empty(&G.graph, nnode, IGRAPH_UNDIRECTED);

    for (int i = 0; i < nnode; ++i) {
        //
        int label = labellist[i];
        G.node_names[i] = label;
    }
    //
    int nbond = psf.bond.size();
    for (int i = 0; i < nbond; ++i) {
        int atm1  = psf.bond[i][0];
        int atm2  = psf.bond[i][1];
        int node1 = find_node_by_label(G, atm1);
        int node2 = find_node_by_label(G, atm2);
        //
        if (node1 != -1 && node2 != -1) {
            igraph_add_edge(&G.graph, node1, node2);
            //
            cout << "Edge added between nodes with labels "\
                 << atm1 << " and " << atm2 << endl;
        }
    }
    //
    return G;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_graph GenerateCommunity(const s_graph& G_bond,\
                          const vector<vector<double>>& coord,\
                          const double& cutoff) {
    //
    s_graph G;
    deep_copy_graph(G, G_bond);
    //
    int nnode = G.node_names.size();

    for (int i = 0; i < nnode; ++i) {
        int atm1 = G.node_names[i];
        for (int j = 0; j < nnode; ++j) {
           int  atm2 = G.node_names[j];
           bool sign = edge_exists(G.graph, i, j);
           //
           if ((i != j) && (sign == false)) {
               vector<double> vec(3, 0.0);
               for (int k = 0; k < 3; ++k) {
                   vec[k] = coord[atm1-1][k] - coord[atm2-1][k];
               }
               double length = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
               if (length <= cutoff) {
                  igraph_add_edge(&G.graph, i, j);
               }
           }
        }
    }

    return G;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_graph GenerateCommunityEdge(const vector<vector<int>>& resinfo,
                              const vector<int>& residinfo,
                              const vector<vector<double>>& coord,
                              const double& cutoff,
                              const vector<vector<double>>& wij) {
    //
    s_graph G;
    igraph_vector_init(&G.weights, 0);
    int nres = resinfo.size();
    igraph_empty(&G.graph, nres, IGRAPH_UNDIRECTED);

    for (int i = 0; i < nres; ++i) {
        if (i != nres-1) {
            igraph_add_edge(&G.graph, i, i+1);
            igraph_vector_push_back(&G.weights, wij[i][i+1]);
        }
        G.node_names.push_back(residinfo[i]);
    }

    //
    for (int i = 0; i < nres; ++i) {
        for (int j = 0; j < nres; ++j) {
            bool esign = edge_exists(G.graph, i, j);
            if (i != j && esign == false) {
                int natm_i = resinfo[i].size();
                int natm_j = resinfo[j].size();
                bool loop_sign = true;
                //
                // For Loop about element of nodes
                for (int k = 0; k < natm_i; ++k) {
                    for (int l = 0; l < natm_j; ++l) {
                        vector<double> vec(3, 0.0);
                        int label_i = resinfo[i][k];
                        int label_j = resinfo[j][l];
                        for (int m = 0; m < 3; ++m) {
                            vec[m] = coord[label_i-1][m] - coord[label_j-1][m];
                        }
                        double length = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                        //
                        if (length < cutoff) {
                            igraph_add_edge(&G.graph, i, j);
                            igraph_vector_push_back(&G.weights, wij[i][j]);
                            loop_sign = false;
                            break;
                        }
                    }
                    if (!loop_sign) break;
                }
            }
        }
    }

    print_graph_summary(G.graph);

    return G;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
s_graph GenerateCommunityLen(const vector<vector<int>>& resinfo, 
                             const vector<int>& residinfo, 
                             const vector<vector<double>>& coord, 
                             const double& cutoff) {
    //
    s_graph G;
    igraph_vector_init(&G.weights, 0);
    int nres = resinfo.size();
    igraph_empty(&G.graph, nres, IGRAPH_UNDIRECTED);
    //
    for (int i = 0; i < nres; ++i) {
        G.node_names.push_back(residinfo[i]);
    }
    
    for (int i = 0; i < nres-1; ++i) {
        for (int j = i+1; j < nres; ++j) {
            int natm_i = resinfo[i].size();
            int natm_j = resinfo[j].size();
            //
            // For loop about elements of nodes
            double len_min = 1000.0;
            bool edge_sign = false;
            for (int k = 0; k < natm_i; ++k) {
                for (int l = 0; l < natm_j; ++l) {
                    vector<double> vec(3, 0.0);
                    int labal_i = resinfo[i][k];
                    int labal_j = resinfo[j][l];
                    for (int m = 0; m < 3; ++m) {
                        //cout << coord[labal_i-1][m] << " " << coord[labal_j-1][m] << endl;
                        vec[m] = coord[labal_i-1][m] - coord[labal_j-1][m];
                    }
                    double length = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                    if ((length < cutoff) && (len_min > length)) {
                        len_min = length;
                        edge_sign = true;
                        //igraph_add_edge(&G.graph, i, j);
                        //igraph_vector_push_back(&G.weights, length);
                    }
                }
            }
            //
            if (edge_sign == true) {
                igraph_add_edge(&G.graph, i, j);
                igraph_vector_push_back(&G.weights, len_min);
            }
        }
    }

    print_graph_summary(G.graph);
    
    return G;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Calc. confusion matrix
vector<vector<int>> confusion_matrix(const vector<int>& memb_refs, 
                                     const vector<int>& memb) {
    map<int, int> refs_map, pred_map;
    int n = memb_refs.size();

    for (int label : memb_refs) refs_map[label]++;
    for (int label : memb) pred_map[label]++;

    vector<vector<int>> matrix(refs_map.size(), vector<int>(pred_map.size(), 0));

    map<int, int> refs_index, pred_index;
    int i = 0;
    for (auto& [label, _] : refs_map) refs_index[label] = i++;
    i = 0;
    for (auto& [label, _] : pred_map) pred_index[label] = i++;

    for (int j = 0; j < n; j++) {
        matrix[refs_index[memb_refs[j]]][pred_index[memb[j]]]++;
    }

    return matrix;
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// NMI calc.
double normalized_mutual_information(const vector<int>& memb_refs, 
                                     const vector<int>& memb) {
    auto matrix = confusion_matrix(memb_refs, memb);
    int n = memb_refs.size();

    vector<int> refs_marginals(matrix.size(), 0), pred_marginals(matrix[0].size(), 0);
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            refs_marginals[i] += matrix[i][j];
            pred_marginals[j] += matrix[i][j];
        }
    }

    double h_refs = 0.0, h_pred = 0.0, mi = 0.0;
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] > 0) {
                double p_ij = static_cast<double>(matrix[i][j]) / n;
                double p_i = static_cast<double>(refs_marginals[i]) / n;
                double p_j = static_cast<double>(pred_marginals[j]) / n;
                mi += p_ij * log(p_ij / (p_i * p_j));
            }
        }
    }

    for (int t : refs_marginals)
        if (t > 0) h_refs -= (static_cast<double>(t) / n) * log(static_cast<double>(t) / n);
    for (int p : pred_marginals)
        if (p > 0) h_pred -= (static_cast<double>(p) / n) * log(static_cast<double>(p) / n);

    return mi / ((h_refs + h_pred) / 2.0);
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// ARI calc.
double adjusted_rand_index(const vector<int>& memb_refs, 
                           const vector<int>& memb) {
    auto matrix = confusion_matrix(memb_refs, memb);
    int n = memb_refs.size();

    vector<int> refs_marginals(matrix.size(), 0), pred_marginals(matrix[0].size(), 0);
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            refs_marginals[i] += matrix[i][j];
            pred_marginals[j] += matrix[i][j];
        }
    }

    double sum_index = 0.0, sum_refs = 0.0, sum_pred = 0.0;
    for (size_t i = 0; i < matrix.size(); i++)
        for (size_t j = 0; j < matrix[i].size(); j++)
            sum_index += matrix[i][j] * (matrix[i][j] - 1) / 2.0;
    for (int t : refs_marginals)
        sum_refs += t * (t - 1) / 2.0;
    for (int p : pred_marginals)
        sum_pred += p * (p - 1) / 2.0;

    double expected_index = sum_refs * sum_pred / (n * (n - 1) / 2.0);
    double max_index = (sum_refs + sum_pred) / 2.0;

    return (sum_index - expected_index) / (max_index - expected_index);
}
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
pair<int, double> ReturnCommNMI(const vector<vector<int>>& refs_memb, 
                                const vector<int>& memb) {
    //
    int    ncomp = refs_memb.size();
    int    comm  = 0;
    double nmi_max = 0.0;

    for (int i = 0; i < ncomp; ++i) {
        double nmi = normalized_mutual_information(refs_memb[i], memb);
        if (nmi > nmi_max) {
            nmi_max = nmi;
            comm = i + 1;
        }
    }

    return {comm, nmi_max};
}
//---------------------------------------------------------------------------//
