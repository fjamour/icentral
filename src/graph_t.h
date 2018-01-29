/* 
 * File:   graph_t.h
 * Author: fuad
 *
 * Created on January 6, 2014, 7:55 AM
 */

#ifndef GRAPH_T_H
#define	GRAPH_T_H

#include <vector>
#include <set>
#include <utility>
#include <string>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <stack>
#include <queue>

#include "types.h"

#define tr1_map_t(value_t) tr1::unordered_map< node_id_t, value_t >

using namespace std;

struct node_t {
    vector<node_id_t>     nbrs_vec;        
};

/* minimum cycle basis for a graph is a set of cycles
 * the representation is a vector of edge lists
 */
struct mcb_t {
    std::vector<cycle_t> cycle_vec;
};

struct graph_t;
struct component_t;


/*
 * graph with nodes having integral indices that shouldn't be
 * between 0 and n-1
 * will be used to store MUCs and subgraphs
 */
struct graph_hash_t {
    typedef tr1::unordered_map<node_id_t, node_t> nodes_map_t;
    nodes_map_t      nodes_map;
    set<edge_t >     edge_set;
    
    //structures to be used by algorithms
    //they are here so they are built once
    tr1_map_t(vector<node_id_t> )       P;
    tr1_map_t(int)                      path_cnt_vec;
    tr1_map_t(int)                      new_path_cnt_vec;
    tr1_map_t(int)                      dist_vec;
    tr1_map_t(double)                   pair_dep_vec;
    tr1_map_t(double)                   new_pair_dep_vec;
    tr1_map_t(int)                      path_cnt_inc_vec;
    tr1_map_t(double)                   sigma_t_map;
    tr1_map_t(double)                   new_sigma_t_map;
    tr1_map_t(bool)                     visited_vec;
    vector<node_id_t>                   S;
    queue<node_id_t>                    Q;
    
    void        insert_node(node_id_t id);
    void        insert_edge(node_id_t src, node_id_t dst);
    void        remove_edge(node_id_t src, node_id_t dst);
    size_t      size();
    void        print_graph(bool with_nodes=true);
    bool        has_edge(node_id_t src, node_id_t dst);
    
    void        fill_graph(graph_t& g);
    
    void        conn_comp_sizes(vector<int>& out_vec);
    
    void        find_conn_comp(vector<graph_hash_t>& out_vec);
    
    void        init_maps();
    
    vector<node_id_t>&  get_nbrs(node_id_t id);
    
    void        find_sssp(node_id_t s, tr1_map_t(int)& out_vec);
    
    template <typename T>
    void        i_fill_map(tr1_map_t(T)& vec, T val);
    
    void        find_pruning_counts_exp(node_id_t src,
                                        node_id_t dst,
                                        int& d0_cnt,
                                        int& d1_cnt,
                                        int& d2_cnt); 
};


/*
 * all functions here assume the given node ids are proper ones from 0 to n-1
 * so, the caller must use 
 */
struct subgraph_t {
    vector<vector<node_id_t> >  nodes_vec;
    set<edge_t >                edge_set;
    
    //maps original node id to id in the the subgraph
    vector<node_id_t>           inout_label_map; // maps internal labels to original
    tr1_map_t(node_id_t)        outin_label_map; // original labels to internal
    
    //structures to be used by algorithms
    //they are here so they are built once
    vector<vector<node_id_t> >       P;
    vector<int>                      path_cnt_vec;
    vector<int>                      new_path_cnt_vec;
    vector<int>                      dist_vec;
    vector<double>                   pair_dep_vec;
    vector<double>                   new_pair_dep_vec;
    vector<int>                      path_cnt_inc_vec;
    vector<double>                   sigma_t_map;
    vector<double>                   new_sigma_t_map;
    vector<bool>                     visited_vec;
    vector<node_id_t>                S;
    queue<node_id_t>                 Q;
    
    void        insert_edge(node_id_t src, node_id_t dst);
    void        remove_edge(node_id_t src, node_id_t dst);
    size_t      size();
    void        print_graph(bool with_nodes=true);
    bool        has_edge(node_id_t src, node_id_t dst);
    
    void        fill_graph(graph_hash_t& g);
    
    void        fill_graph(graph_t& g);
    
    void        conn_comp_sizes(vector<int>& out_vec);
    
    void        init_maps();
    
    vector<node_id_t>&  get_nbrs(node_id_t id);
    
    void        find_sssp(node_id_t s, vector<int>& out_vec);
    
    template <typename T>
    void        i_fill_map(vector<T>& vec, T val){
        vec.resize(nodes_vec.size());
        fill(vec.begin(), vec.end(), val);
    }
    
    void        find_pruning_counts_exp(node_id_t src,
                                        node_id_t dst,
                                        int& d0_cnt,
                                        int& d1_cnt,
                                        int& d2_cnt); 
};

/* not sure yet TODO
 */
struct muc_t {
    typedef tr1::unordered_map<node_id_t, vector<node_id_t> >   conn_vert_map_t;
    typedef tr1::unordered_map<node_id_t, graph_hash_t>         subgraph_map_t;
    
    conn_vert_map_t     conn_vertex_map;
    subgraph_map_t      subgraph_map;
    graph_hash_t        muc_subgraph;
    muc_id_t            id;
    bool                valid;//flag to tell if the muc was deleted or not
    
    //these are used to facilitate fast computation where iterations are done in
    //a fast graph, not hash based one
    subgraph_t          muc_fast_subgraph;
    conn_vert_map_t     tmp_conn_vertex_map;
    subgraph_map_t      tmp_subgraph_map;
    
    muc_t();
    
    void        insert_conn_vertex(node_id_t, node_id_t nbr);
    void        compute_bc(tr1::unordered_map<node_id_t, double>& bc_map,
                           int max_iter = -1);
    
    
    //incremental BFS stuff below
    void        compute_bc_inc(tr1_map_t(double)& bc_map,
                        node_id_t src, 
                        node_id_t dst,
                        int max_iter = -1);
    
    
    void        i_iteration(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_vec);
    
    //|d_src-d_dst| = 1 (the easy case) 
    void        i_iteration_1(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_vec);
    //|d_src-d_dst| >= 2 (the difficult case)
    void        i_iteration_2(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_vec);
};



struct bcc_scratch_t {
    //maps articulation points to sizes of subgraphs connected to the bcc
    //through them
    typedef tr1::unordered_map<node_id_t, vector<int> >   art_pt_map_t;
    
    art_pt_map_t        art_pt_map;
    graph_hash_t        bcc_subgraph;
    subgraph_t          bcc_fast_subgraph;
   
    //This compute_bc DOES NOT return correct bc values for articulation points
    void        compute_bc(tr1::unordered_map<node_id_t, double>& bc_map, int max_iter = -1);
    void        bc_iter(node_id_t s, tr1::unordered_map<node_id_t, double>& bc_map);
    void        print();
};

struct bcc_delta_t {
    //maps articulation points to sizes of subgraphs connected to the bcc
    //through them
    typedef tr1_map_t(vector<int>)   art_pt_map_t;
    
    art_pt_map_t        art_pt_map;
    graph_hash_t        bcc_subgraph;
    subgraph_t          bcc_fast_subgraph;
   
    void        compute_bc(tr1_map_t(double)& bc_map,
                        node_id_t src, 
                        node_id_t dst);
    void        compute_bc_exp(tr1_map_t(double)& bc_map,
                        node_id_t src, 
                        node_id_t dst,
                        bcc_stat_t& bcc_stat,
                        int max_iter_d1 = -1,
                        int max_iter_d2 = -1);
    
    void        compute_bc_maxiter_exp(tr1_map_t(double)& bc_map,
                        node_id_t src,
                        node_id_t dst,
                        bcc_stat_t& bcc_stat,
                        int max_iter_d1,
                        int max_iter_d2);
    void        print();
    
    void        i_iteration(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_vec);
    
    //|d_src-d_dst| = 1 (the easy case) 
    void        i_iteration_1(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_vec);
    //|d_src-d_dst| >= 2 (the difficult case)
    void        i_iteration_2(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_vec);
    
    void        dbg_iteration(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_vec);
    
};

/*
 * Simple undirected unweighted graph data structure
 * no checks whatsoever
 * should call init_size(..) first then insert_edge(..) to populate
 * 
 * IMP: nodes have indexes from 0 to n-1
 * IMP: graph is assumed to be connected
 */
struct graph_t {
    vector<node_t>      nodes_vec;
    vector<muc_id_t>    node_to_muc_vec;
    set<edge_t >        edge_set;
    vector<muc_t>       muc_vec;
    mcb_t               mcb;
    
    vector<bool>        tmp_bool_vec;
    
    bool                bc_computed;
    vector<double>      bc_vec;
    tr1_map_t(double)   bc_map;
    
    string              graph_name;

    void                insert_edge(node_id_t src, node_id_t dst);
    void                remove_edge(node_id_t src, node_id_t dst);
    void                init_size(size_t num_nodes);
    
    size_t              size();
    void                read_graph(string path);
    vector<node_id_t>&  get_nbrs(node_id_t id);
    void                get_shortest_path(
                                node_id_t src,
                                node_id_t dst,
                                vector<node_id_t>& out_path);
    bool                has_edge(const edge_t& e);
    
    void                print_edgelist();
    void                print_mcbs();
    void                print_mucs(bool with_nodes=true);
    void                print_node_to_muc();
    
    void                find_mucs();
    void                find_mucs_fast();
    void                find_muc_mcb();
    void                find_conn_verts();
    void                merge_muc_cycle(muc_t& muc, cycle_t& cycle);
    void                merge_muc_muc(muc_t& muc1, muc_t& muc2);
    void                find_muc_subgraphs(muc_id_t muc_id);
    size_t              get_num_mucs();
    
    void                insert_edge_update_muc(node_id_t src, node_id_t dst);
    
    
    
    //bcc stuf
    void                find_art_points(vector<node_id_t>& out_vec);
    void                find_bridge_edges(vector<edge_t>& out_vet);
    void                find_bicon_comp(vector<graph_hash_t>& out_vec);
    void                find_edge_bcc_subgraph(graph_hash_t& bcc,
                                node_id_t src,
                                node_id_t dst);
    
    void                find_edge_bcc(bcc_scratch_t& bcc,
                                node_id_t src,
                                node_id_t dst);
    
    void                find_edge_bcc(bcc_delta_t& bcc,
                                node_id_t src,
                                node_id_t dst);
    
    void                find_sssp(node_id_t s, vector<int>& out_vec);
    
    
    void                tmp_fun();
    
    void                i_init_internals();
    
    void                init_bc();
    
    void    find_edge_bcc(
                component_t&    comp,
                edge_t          e
                );
    
    //incremental bc with bcc stuff
    void                insert_edge_update_bc_experimental(node_id_t src,
                                                           node_id_t dst,
                                                           bcc_stat_t& bcc_stat,
                                                           int max_iter_d1 = -1,
                                                           int max_iter_d2 = -1);
   
    //approximates the runtime of Brandes iteration on a bcc subgraph
    //by doing @num_iter iteration and taking the average
    //if @num_iter is -1 the number of iteration is just going to be the
    //number of nodes in the graph
    void                approx_bcc_iter_tm(node_id_t src,
                                        node_id_t dst,
                                        double& tm,
                                        int num_iter = -1);
private:
    
    
    void                i_art_point_dfs_visitor(node_id_t u,
                                    vector<int>& color_vec,
                                    vector<int>& low_vec,
                                    vector<int>& d_vec,
                                    vector<node_id_t>& pred_vec,
                                    int& time,
                                    vector<node_id_t>& art_point_vec);
    void                i_bcc_dfs_visitor(node_id_t u,
                                    vector<int>& color_vec,
                                    vector<int>& low_vec,
                                    vector<int>& d_vec,
                                    vector<node_id_t>& pred_vec,
                                    stack<edge_t>& edge_stack,
                                    int& time,
                                    vector<graph_hash_t>& bcc_vec);
    void                i_edge_bcc_dfs_visitor(node_id_t u,
                                    vector<int>& color_vec,
                                    vector<int>& low_vec,
                                    vector<int>& d_vec,
                                    vector<node_id_t>& pred_vec,
                                    stack<edge_t>& edge_stack,
                                    int& time);
};





#endif	/* GRAPH_T_H */

