/* 
 * File:   unit_tests.h
 * Author: fuad
 *
 * Created on December 22, 2014, 11:21 AM
 */

#ifndef UNIT_TESTS_H
#define	UNIT_TESTS_H

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <list>
#include <algorithm>
#include <cmath>


#include "graph_t.h"
#include "bc.h"
#include "experiments.h"
#include "utility.h"
#include "types.h"

#include <mpi.h>

void fill_test_graph(graph_t& graph)
{
    graph.init_size(20);
    graph.insert_edge(0,1);
    graph.insert_edge(1,2);
    graph.insert_edge(2,3);
    graph.insert_edge(3,0);
    graph.insert_edge(0,4);
    graph.insert_edge(4,5);
    graph.insert_edge(5,6);
    graph.insert_edge(6,7);
    graph.insert_edge(7,8);
    graph.insert_edge(8,4);
    graph.insert_edge(5,9);
    graph.insert_edge(9,8);
    graph.insert_edge(3,10);
    graph.insert_edge(10,11);
    graph.insert_edge(11,12);
    graph.insert_edge(12,13);
    graph.insert_edge(13,14);
    graph.insert_edge(14,15);
    graph.insert_edge(15,10);
    graph.insert_edge(2,16);
    graph.insert_edge(16,17);
    graph.insert_edge(17,18);
    graph.insert_edge(2,19);
}

bool test__fast_brandes_BC()
{
    bool passed = false;
    printf("Testing [brandes_BC] ...\t");
    graph_t g;
    g.read_graph(string("Erdos02.lcc.net"));
    //fill_test_graph(g);
    vector<double> BC_vec;
    fast_brandes_BC(g, BC_vec);
    vector<double> ref_BC_vec = brandes_bc(g);
    passed = true;
    for(node_id_t v = 0; v < BC_vec.size(); ++v) {
        if((BC_vec[v] - ref_BC_vec[v]) > EPS) {
            passed = false;
        }
        //printf("[%f]\t\t[%f]\n", BC_vec[v], ref_BC_vec[v]);
    }
    if(passed) {
        printf("PASSED\n");
    } else {
        printf("FAILED\n");
    }
    
    return passed;
}

bool test__iCentral()
{
    bool passed = false;
    printf("Testing [iCentral] ...\t");
    
    edge_t e = make_pair(3, 18);
    graph_t g;
    g.read_graph(string("Erdos02.lcc.net"));
    //fill_test_graph(g);
    vector<double> BC_vec;
    fast_brandes_BC(g, BC_vec);
    
    component_t comp;
    comp.comp_type = GRAPH;
    comp.subgraph.fill_graph(g);
    vector<double> dBC_vec;
    iCentral(dBC_vec, comp, e);
    for(int i = 0; i < dBC_vec.size(); ++i) {
        BC_vec[i] += dBC_vec[i];
    }
    
    g.insert_edge(e.first, e.second);
    vector<double> ref_BC_vec = brandes_bc(g);
    passed = true;
    for(node_id_t v = 0; v < BC_vec.size(); ++v) {
        if(abs((BC_vec[v] - ref_BC_vec[v])) > EPS) {
            passed = false;
        }
        //printf("[%f]\t\t[%f]\n", BC_vec[v], ref_BC_vec[v]);
    }
    if(passed) {
        printf("PASSED\n");
    } else {
        printf("FAILED\n");
    }
    
    return passed;
}

bool test__Update_BC()
{
    bool passed = false;
    printf("Testing [Update_BC] ...\t");
    vector<double> BC_vec;
    graph_t g;
    edge_t e = make_pair(100, 3455);
    
    g.read_graph(string("Erdos02.lcc.net"));
    //fill_test_graph(g);
    
    e = *g.edge_set.begin();
    g.remove_edge(e.first, e.second);
    //printf("(%d, %d)\n", e.first, e.second);
       
    fast_brandes_BC(g, BC_vec);
    Update_BC(BC_vec, g, BCC, e);
    
    g.insert_edge(e.first, e.second);
    vector<double> ref_BC_vec = brandes_bc(g);
    passed = true;
    for(node_id_t v = 0; v < BC_vec.size(); ++v) {
        if(abs((BC_vec[v] - ref_BC_vec[v])) > EPS) {
            passed = false;
            printf("(%d) --\t[%f]\t\t[%f]\n", v, BC_vec[v], ref_BC_vec[v]);
        }
        //printf("(%d) --\t[%f]\t\t[%f]\n", v, BC_vec[v], ref_BC_vec[v]);
    }
    if(passed) {
        printf("PASSED\n");
    } else {
        printf("FAILED\n");
    }
    return passed;
}


/*
 * Actual timing details are done here
 */
void timing__Update_BC_graph(
            graph_t&        graph,
            vector<edge_t>& edge_vec,
            int             num_sources,
            comp_type_t     algo_flag,
            bool            do_brandes = true,
            double          brandes_time = 0.0,
            int             num_threads = 1,
            bool            del_edge = false,//in case the edges are already in the graph
            operation_t     op = INSERTION
        )
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status    status;
    
    timer           tm;
    vector<double>  BC_vec;
    vector<double>  tm_vec;
    vector<double>  speedup_vec;
    
    BC_vec.resize(graph.size());
    if(do_brandes) {
        tm.start();
        //fast_brandes_BC(graph, BC_vec);
        BC_vec = brandes_bc(graph);
        tm.stop();
        brandes_time = tm.interval();
    }
    if(rank == 0) {
        printf("Graph[%s]  V[%d]  E[%d]  Brandes_tm[%.2f]\n",
                graph.graph_name.c_str(),
                graph.size(),
                graph.edge_set.size(),
                brandes_time);
    }
    for(int i = 0; i < edge_vec.size(); ++i) {
        edge_t e = edge_vec[i];
        if(del_edge) graph.remove_edge(e.first, e.second);
        tm.start();
        Update_BC(BC_vec, graph, algo_flag, e, num_threads, op);
        tm.stop();
        double e_time = tm.interval();
        tm_vec.push_back(e_time);
        double e_speedup = brandes_time/e_time;
        speedup_vec.push_back(e_speedup);
        
        if(rank == 0) {
            printf("e(%-6d,%-6d)  tm[%.2f]  sup[%.2f]\n",
                    e.first,
                    e.second,
                    e_time,
                    e_speedup);
        }
        if(del_edge) graph.insert_edge(e.first, e.second);
        //synchronization barrier so no one starts next edge before others
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    double tm_mean, tm_stddev, tm_median, tm_min, tm_max;
    double su_mean, su_stddev, su_median, su_min, su_max;
    simple_stats(tm_vec, tm_mean, tm_stddev, tm_median, tm_min, tm_max);
    simple_stats(speedup_vec, su_mean, su_stddev, su_median, su_min, su_max);
    
    if(rank == 0)
        printf("Avg.tm[%.2f]  Avg.sup[%.2f]\n\n", tm_mean, su_mean);
}


/*
 * Does timing test for graphs in @path_vec
 */
void timing__Update_BC(
            vector<string>  path_vec,           //paths of the graphs
            int             num_edges,          //number of random edges
            int             rand_seed,          //seed for srand(...)
            int             num_sources,        //
            comp_type_t     algo_flag = BCC,    //what algorithms to evaluate
            int             num_thread = 1,
            operation_t     op = INSERTION
        )
{
    srand(rand_seed);
    vector<vector<edge_t> > edge_vec2;
    
    //generate random edges
    for(int p = 0; p < path_vec.size(); ++p) {
        string graph_path = path_vec[p];
        graph_t graph;
        graph.read_graph(graph_path);
        graph.graph_name = extract_graph_name(graph_path);
        vector<edge_t> edge_vec;
        gen_rand_edges(num_edges, graph, edge_vec);
        edge_vec2.push_back(edge_vec);
    }
    
    for(int p = 0; p < path_vec.size(); ++p) {
        string graph_path = path_vec[p];
        graph_t graph;
        graph.read_graph(graph_path);
        graph.graph_name = extract_graph_name(graph_path);
        timing__Update_BC_graph(graph,
                edge_vec2[p], 
                num_sources, 
                algo_flag, 
                false, 
                1.0,
                num_thread,
                op
                );
    }
}



/*************************************************
 * EXPERIMENTAL
 *************************************************/

void timing__Update_BC_mem_graph(
            graph_t&        graph,
            vector<edge_t>& edge_vec,
            int             num_sources,
            int             algo_flag,
            bool            do_brandes = true,
            double          brandes_time = 0.0
        )
{
    timer           tm;
    vector<double>  BC_vec;
    vector<double>  tm_vec;
    vector<double>  speedup_vec;
    
    
    if(do_brandes) {
        tm.start();
        //fast_brandes_BC(graph, BC_vec);
        BC_vec = brandes_bc(graph);
        tm.stop();
        brandes_time = tm.interval();
    }
    
    
    
    vector<iter_info_t> iter_info_vec;
    component_t comp;
    comp.comp_type = GRAPH;
    comp.subgraph.fill_graph(graph);
    for (node_id_t s = 0; s < graph.size(); ++s) {
        iter_info_t iter_info;
        iter_info.init_all(graph.size());
        BBFS(iter_info, comp, s);
        iter_info_vec.push_back(iter_info);
    }
           
    printf("Graph[%s]  V[%d]  E[%d]  Brandes_tm[%.2f]\n",
            graph.graph_name.c_str(),
            graph.size(),
            graph.edge_set.size(),
            brandes_time);
    for(int i = 0; i < edge_vec.size(); ++i) {
        edge_t e = edge_vec[i];
        tm.start();
        Update_BC_mem(BC_vec, graph, GRAPH, e, iter_info_vec);
        tm.stop();
        double e_time = tm.interval();
        tm_vec.push_back(e_time);
        double e_speedup = brandes_time/e_time;
        speedup_vec.push_back(e_speedup);
        
        printf("e(%-6d,%-6d)  tm[%.2f]  sup[%.2f]\n",
                e.first,
                e.second,
                e_time,
                e_speedup);
    }
    
    double tm_mean, tm_stddev, tm_median, tm_min, tm_max;
    double su_mean, su_stddev, su_median, su_min, su_max;
    simple_stats(tm_vec, tm_mean, tm_stddev, tm_median, tm_min, tm_max);
    simple_stats(speedup_vec, su_mean, su_stddev, su_median, su_min, su_max);
    
    printf("Avg.tm[%.2f]  Avg.sup[%.2f]\n\n", tm_mean, su_mean);
}

void timing__Update_BC_mem(
            vector<string>  path_vec,       //paths of the graphs
            int             num_edges,      //number of random edges
            int             rand_seed,      //seed for srand(...)
            int             num_sources,    //
            int             algo_flag       //what algorithms to evaluate
        )
{
    vector<vector<edge_t> > edge_vec2;
    
    //generate random edges
    for(int p = 0; p < path_vec.size(); ++p) {
        string graph_path = path_vec[p];
        graph_t graph;
        graph.read_graph(graph_path);
        graph.graph_name = extract_graph_name(graph_path);
        vector<edge_t> edge_vec;
        gen_rand_edges(num_edges, graph, edge_vec);
        edge_vec2.push_back(edge_vec);
    }
    
    for(int p = 0; p < path_vec.size(); ++p) {
        string graph_path = path_vec[p];
        graph_t graph;
        graph.read_graph(graph_path);
        graph.graph_name = extract_graph_name(graph_path);
        timing__Update_BC_mem_graph(graph, edge_vec2[p], -1, 1);
    }
}


bool test__Update_BC_mem()
{
    bool passed = false;
    printf("Testing [Update_BC_mem] ...\t");
    vector<double> BC_vec;
    graph_t g;
    
    
    edge_t e = make_pair(0, 2);
    //fill_test_graph(g);
    g.read_graph(string("Erdos02.lcc.net"));
    e = *(g.edge_set.begin());
    g.remove_edge(e.first, e.second);
    
    printf("(%d, %d)\n", e.first, e.second);
    
    
    vector<iter_info_t> iter_info_vec;
    component_t comp;
    comp.comp_type = GRAPH;
    comp.subgraph.fill_graph(g);
    vector<double> dummy_BC_vec;
    dummy_BC_vec.resize(g.size());
    for (node_id_t s = 0; s < g.size(); ++s) {
        iter_info_t iter_info;
        iter_info.init_all(g.size());
        BBFS(iter_info, comp, s);
        RBFS(dummy_BC_vec, comp, s, iter_info, false, false);
        iter_info_vec.push_back(iter_info);
    }
    
    fast_brandes_BC(g, BC_vec);
    Update_BC_mem(BC_vec, g, GRAPH, e, iter_info_vec);
    
    g.insert_edge(e.first, e.second);
    vector<double> ref_BC_vec = brandes_bc(g);
    passed = true;
    for(node_id_t v = 0; v < BC_vec.size(); ++v) {
        if(abs((BC_vec[v] - ref_BC_vec[v])) > EPS) {
            passed = false;
            printf("(%d) --\t[%f]\t\t[%f]\n", v, BC_vec[v], ref_BC_vec[v]);
        }
        //printf("(%d) --\t[%f]\t\t[%f]\n", v, BC_vec[v], ref_BC_vec[v]);
    }
    if(passed) {
        printf("PASSED\n");
    } else {
        printf("FAILED\n");
    }
    return passed;
}

#endif	/* UNIT_TESTS_H */

