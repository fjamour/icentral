/* 
 * File:   experiments.h
 * Author: fuad
 *
 * Created on March 2, 2014, 10:29 AM
 */

#ifndef EXPERIMENTS_H
#define	EXPERIMENTS_H

#include "graph_t.h"

struct edge_stat_t {
    double      edge_ins_time;
    double      muc_bc_update_time;
    int         muc_num_nodes;
    int         muc_num_edges;
};

void gen_rand_edges(int num_edges,
        graph_t& graph,
        vector<edge_t>& out_vec);

void gen_rand_edges_deletions(int num_edges,
        graph_t& graph,
        vector<edge_t>& out_vec);

void insertion_test_qube(graph_t& graph, int num_edges_to_insert=1);

void incremental_brandes_test(graph_t& graph,
        int num_edges_to_insert=1,
        vector<edge_t> rand_edge_vec = vector<edge_t>());

void insertion_test_qube_hash(graph_t& graph,
        int num_edges_to_insert,
        vector<edge_t> rand_edge_vec = vector<edge_t>(),
        bool use_incremental_algo = false,
        int max_iter = -1);

void insertion_test_fuad_hash(graph_t& graph, int num_edges_to_insert);

void insertion_test_fuad_max_iter(graph_t& graph,
        int num_edges_to_insert,
        int num_iter,
        double time_to_run_all_iter = -1,
        vector<edge_t> rand_edge_vec = vector<edge_t>());

void compare_brandes_hash_vs_vector(graph_t& graph, int num_iter);

void speedup_info_lbcc(graph_t& graph, int num_edges_to_insert);

void speedup_info(graph_t& graph, int num_edges_to_insert);

void count_bcc(graph_t& graph, int num_edges);

void qube_ideal_speedup(graph_t& graph,
        int num_edges_to_insert,
        vector<edge_t> rand_edge_vec = vector<edge_t>());

void lbcc_stat(graph_t& graph);


//paper experiments
double  exp_brandes_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double time_to_run_all_iter = -1);
void    exp_inc_brandes_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double brandes_time);
void    exp_qube_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double brandes_time);
void    exp_inc_qube_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double brandes_time);
void    exp_fuad_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double brandes_time);



#endif	/* EXPERIMENTS_H */

