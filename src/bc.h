/* 
 * File:   bc.h
 * Author: fuad
 *
 * Created on January 6, 2014, 9:27 AM
 */

#ifndef BC_H
#define	BC_H

#include "graph_t.h"

#include "iter_info_t.h"
#include "component_t.h"

using namespace std;

enum operation_t {INSERTION, DELETION};

vector<double> naive_betweenness_centrality(graph_t& graph_in);

vector<double> brandes_betweenness_centrality(graph_t& graph_in);

vector<double> brandes_bc(graph_t& graph_in, int max_iter = -1);

tr1_map_t(double) brandes_bc_hash_out(graph_t& graph_in, int max_iter = -1);

tr1_map_t(double) brandes_bc_hash(graph_hash_t& graph_in, int max_iter = -1);

//increments bc values in @bc_vec, with the pair dependency of @s
void brandes_iter(graph_t& graph_in,
        node_id_t s,
        vector<double>& bc_vec);

//increments bc values in @bc_vec, with the pair dependency of @s
void brandes_iter_hash_out(graph_t& graph_in,
        node_id_t s,
        tr1_map_t(double)& bc_vec);

//increments bc values in @bc_vec, with the pair dependency of @s
void brandes_iter(graph_hash_t& graph_in,
        node_id_t s,
        tr1_map_t(double)& bc_vec);


void incremental_brandes(graph_t& graph_in,
        node_id_t src,
        node_id_t dst,
        vector<double>& bc_vec);

void incremental_brandes_experimental(graph_t& graph_in,
        node_id_t src,
        node_id_t dst,
        vector<double>& bc_vec,
        vector<double>& time_vec,
        vector<int>&    cnt_vec);


//increments bc values in @bc_vec with the change of pair dependency of @s
//after adding edge (src, dst)
//edge (src, dst) is assumed not to be in the graph
void brandes_delta_iter(graph_t& graph_in,
        node_id_t s,
        node_id_t src,
        node_id_t dst,
        int d_src,
        int d_dst,
        vector<double>& bc_vec);




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
/*
 * New functions
 */

/*
 * Computes betweenness centrality using Brandes algorithm
 * Will construct @BC_vec
 */
void fast_brandes_BC(
        graph_t& graph,
        vector<double>& BC_vec
        );


/*
 * Updates the BC values in @BC_vec in place for nodes in @graph using
 * @comp_type decomposition
 * TODO: should handle all kinds of decompositions (graph/BCC/MUC)
 * IMP: assumes @BC_vec is the right size and has BC values
 */
void Update_BC(
        vector<double>&  BC_vec,
        graph_t&         graph,
        comp_type_t      comp_type,
        edge_t           e,
        int              num_threads = 1,
        operation_t      op = INSERTION
        );

/*
 * Computes increments to BC in in @comp.subgraph after edge @e is inserted
 * TODO: should handle any kind of component (graph/BCC/MUC)
 * 
 */
void iCentral(
        vector<double>& dBC_vec,     //values will be updated in place
        component_t&    comp,
        edge_t          e,
        operation_t     op = INSERTION
        );

/*
 * Computes the increments/decrements to BC of a subgraph in @comp
 * This function deals with nodes indexed from 0 to N-1 in the passed subgraph
 * and knows nothing about the original graph, the caller must add the
 * deltas to the BC vector of the original graph
 */
void iCentral_iter(
        vector<double>& dBC_vec,    // delta BC of vertices
        component_t&    comp,       // component could be BCC, MUC, or just a graph
        node_id_t       s,          // source of the iteration
        edge_t          e,          // inserted edge
        iter_info_t&    iter_info,   //TODO to be used later?
        int             dd = -1,
        bool            use_d_1 = true,
        operation_t     op = INSERTION
        );

/*
 * One iteration of Brandes algorithm from source @s
 */
void brandes_iter(
        vector<double>& BC_vec,     // BC of vertices
        component_t&    comp,       // component could be BCC, MUC, or just a graph
        node_id_t       s,          // source of the iteration
        iter_info_t&    iter_info  //TODO to be used later?
        );

/*
 * If dd=1 will compute both old values and new values to eliminate need for
 * PartialBFS later
 * IMP: e.first is assumed to be closer to the source than e.second
 */
void BBFS(
        iter_info_t&    iter_info,  // iteration info to be computed
        component_t&    comp,       // component
        node_id_t       s          // source of the iteration
        );

/*
 *
 */
void RBFS(
        vector<double>& dBC_vec,    // delta BC of vertices
        component_t&    comp,       // component could be BCC, MUC, or just a graph
        node_id_t       s,          // source of the iteration
        iter_info_t&    iter_info,  // 
        bool            add,
        bool            sub
        );
/*
 *
 */
void partial_BBFS_addition(
        iter_info_t&    iter_info,  // iteration info to be computed
        component_t&    comp,       // component
        node_id_t       s,          // source of the iteration
        edge_t          e           // inserted edge
        );

void partial_BBFS_deletion(
        iter_info_t&    iter_info,  // iteration info to be computed
        component_t&    comp,       // component
        node_id_t       s,          // source of the iteration
        edge_t          e           // deleted edge
        );

void BBFS_RBFS_d1(
        vector<double>& dBC_vec,    // delta BC of vertices
        iter_info_t&    iter_info,  // iteration info to be computed
        component_t&    comp,       // component
        node_id_t       s,          // source of the iteration
        edge_t          e
        );


void parallel_iCentral(
        vector<double>& dBC_vec,
        component_t&    comp,
        edge_t          e,
        int             num_threads,
        operation_t     op = INSERTION
        );


void iCentral_block(
        vector<double>*     dBC_vec,
        component_t*        comp,
        edge_t*             e,
        vector<node_id_t>*  source_vec,
        operation_t*        op
        );



/*************************************************
 * EXPERIMENTAL
 *************************************************/
void Update_BC_mem(
        vector<double>&         BC_vec,
        graph_t&                graph,
        comp_type_t             comp_type,
        edge_t                  e,
        vector<iter_info_t>&    iter_info_vec
        );

void iCentral_mem(
        vector<double>&         dBC_vec,
        component_t&            comp,
        edge_t                  e,
        vector<iter_info_t>&    iter_info_vec
        );
/*
 * Computes the increments/decrements to BC of a subgraph in @comp
 * This function deals with nodes indexed from 0 to N-1 in the passed subgraph
 * and knows nothing about the original graph, the caller must add the
 * deltas to the BC vector of the original graph
 * Assumes @iter_info has shortest path info, so no BBFS is done
 * IMP: for now works for GRAPH comp only
 */
void iCentral_iter_mem(
        vector<double>& dBC_vec,    // delta BC of vertices
        component_t&    comp,       // component could be BCC, MUC, or just a graph
        node_id_t       s,          // source of the iteration
        edge_t          e,          // inserted edge
        iter_info_t&    iter_info   //TODO to be used later?
        );




/*
 * Parallel Brandes functions
 * Ziyad territory 
 */
void parallel_brandes(
        graph_t&        graph,
        vector<double>& BC_vec
        );

void brandes_block(
        vector<double>*     dBC_vec,
        component_t*        comp,
        vector<node_id_t>*  source_vec
        );

#endif	/* BC_H */
