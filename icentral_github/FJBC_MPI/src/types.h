/* 
 * File:   types.h
 * Author: fuad
 *
 * Created on January 6, 2014, 8:30 AM
 */

#ifndef TYPES_H
#define	TYPES_H

#include <utility>

#define INF 100000000
#define EPS 0.0001

typedef unsigned int                        node_id_t;
typedef std::pair<node_id_t, node_id_t>     edge_t;

typedef unsigned int                        muc_id_t;//TODO remove
typedef std::vector<edge_t>                 cycle_t;//TODO remove

struct bcc_stat_t {//TODO remove
    int         bcc_num_nodes;//size of the bcc this edge insertion
    int         bcc_num_edges;
    double      avg_iter_tm_brandes;//time it takes to do Brandes straight forward iteration 
    int         num_d0_iter;
    int         num_d1_iter;
    int         num_d2_iter;
    double      tot_d0_tm;
    double      tot_d1_tm;
    double      tot_d2_tm;
    double      bcc_find_time;
    double      bc_update_time;
    double      sssp_tm;
    
    int         tot_d0_iter;
    int         tot_d1_iter;
    int         tot_d2_iter;
};


#endif	/* TYPES_H */

