#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <stack>
#include <queue>
#include <numeric>
#include <cmath>

#include "bc.h"
#include "utility.h"

#include <thread>
#include <mpi.h>

using namespace std;

void fast_brandes_BC(
        graph_t&         graph,
        vector<double>& BC_vec
        )
{
    //1. make a component (fill a subgraph)
    //2. loop through all sources
    BC_vec.resize(graph.size());
    fill(BC_vec.begin(), BC_vec.end(), 0.0);
    
    component_t comp;
    comp.comp_type = GRAPH;
    comp.subgraph.fill_graph(graph);
    
    iter_info_t iter_info;
    for(node_id_t s = 0; s < comp.subgraph.size(); ++s) {
        brandes_iter(BC_vec, comp, s, iter_info);
    }
}


/*
 * For edge additions, the edge must not be in the graph
 * For edge deletions, the edge must not be a bridge
 */
void Update_BC(
        vector<double>&  BC_vec,
        graph_t&         graph,
        comp_type_t      comp_type,
        edge_t           e,
        int              num_threads,
        operation_t      op
        )
{
    if(comp_type == GRAPH) {
        component_t comp;
        comp.comp_type = GRAPH;
        comp.subgraph.fill_graph(graph);
        vector<double> dBC_vec;
        if(num_threads == 1)
            iCentral(dBC_vec, comp, e, op);
        else
            parallel_iCentral(dBC_vec, comp, e, num_threads, op);
        for(int i = 0; i < dBC_vec.size(); ++i) {
            BC_vec[i] += dBC_vec[i];
        }
    } else if(comp_type == BCC) {
        component_t comp;
        comp.comp_type = BCC;
        if(op == INSERTION) {
            //IMP: assumes @e is not in @graph, @e will not be in @comp
            graph.find_edge_bcc(comp, e);
        } else if(op == DELETION) {
            
            //IMP: assumes @e is not in @graph, @e will not be in @comp
            //     (remove @e from @graph)
            graph.remove_edge(e.first, e.second);
            graph.find_edge_bcc(comp, e);
            graph.insert_edge(e.first, e.second);
        }
            
        //DEBUG
        //printf("BCC # vertices [%d], # edges [%d] \n", comp.subgraph.size(), comp.subgraph.edge_set.size());
        
        //map @e to the new edge
        e.first  = comp.subgraph.outin_label_map[e.first];
        e.second = comp.subgraph.outin_label_map[e.second];
        if(op == DELETION) {
            comp.subgraph.insert_edge(e.first, e.second);
        }
        vector<double> dBC_vec;
        //comp.print();
        if(num_threads == 1)
            iCentral(dBC_vec, comp, e, op);
        else
            parallel_iCentral(dBC_vec, comp, e, num_threads, op);
        for(int i = 0; i < dBC_vec.size(); ++i) {
            node_id_t actual_node_id = comp.subgraph.inout_label_map[i];
            BC_vec[actual_node_id] += dBC_vec[i];
        }
    } else if(comp_type == MUC) {
        //TO BE IMPLEMENTED
    }
   
    
}

void iCentral(
        vector<double>& dBC_vec,
        component_t&    comp,
        edge_t          e,
        operation_t     op
        )
{
    fill_vec<double>(dBC_vec, comp.subgraph.size(), 0.0);
    vector<int> d_src_vec, d_dst_vec;
    comp.subgraph.find_sssp(e.first, d_src_vec);
    comp.subgraph.find_sssp(e.second, d_dst_vec);
    iter_info_t iter_info;
    for(node_id_t s = 0; s < comp.subgraph.size(); ++s) {
        if(d_src_vec[s] != d_dst_vec[s]) {
            int dd = d_src_vec[s] - d_dst_vec[s]; //dd=d_v1-d_v2
            iCentral_iter(dBC_vec, comp, s, e, iter_info, dd, false ,op);
        }
    }
}


void parallel_iCentral(
        vector<double>& dBC_vec,
        component_t&    comp,
        edge_t          e,
        int             num_threads,
        operation_t     op
        )
{
    fill_vec<double>(dBC_vec, comp.subgraph.size(), 0.0);
    vector<int> d_src_vec, d_dst_vec;
    comp.subgraph.find_sssp(e.first, d_src_vec);
    comp.subgraph.find_sssp(e.second, d_dst_vec);
    vector<node_id_t> all_sources_vec;
    for(node_id_t s = 0; s < comp.subgraph.size(); ++s) {
        if(d_src_vec[s] != d_dst_vec[s]) {
            all_sources_vec.push_back(s);
        }
    }
    //MPI shit goes here:
    //all_sources_vec at each machine must have only its share of the sources
    //then continue normally, note that each process will finish its shit
    //and have it's contribution in its dBC_vec
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status    status;
    if(all_sources_vec.size() < size)
        return;
    vector<node_id_t> machine_source_vec;
    int num_s_per_machine = all_sources_vec.size()/size;
    int start, end;
    start   = rank * num_s_per_machine;
    end     = start + num_s_per_machine;
    if(rank == size - 1)
        end = all_sources_vec.size();
    for(int i = start; i < end; ++i) {
        machine_source_vec.push_back(all_sources_vec[i]);
    }
    all_sources_vec = machine_source_vec;
    //printf("RANK[%d] -- num sources [%d]\n", rank, all_sources_vec.size());
    //////////////
    vector<vector<node_id_t> > thread_source_vec;
    thread_source_vec.resize(num_threads);
    if(all_sources_vec.size() < num_threads)
        return;
    int num_s_per_thread = all_sources_vec.size()/num_threads;
    int t = -1;
    for(int i = 0; i < all_sources_vec.size(); ++i) {
        if(i%num_s_per_thread == 0 && t < num_threads-1)
            t++;
        thread_source_vec[t].push_back(all_sources_vec[i]);
    }
    
    vector<std::thread> thread_vec;
    thread_vec.resize(num_threads);
    //start the threads
    vector<vector<double> > dBC_vec_vec;
    dBC_vec_vec.resize(num_threads);
    
    for(int t = 0; t < num_threads; ++t) {
        thread_vec[t] = std::thread(iCentral_block, &(dBC_vec_vec[t]), &comp, &e, &(thread_source_vec[t]), &op);
    }
    //wait for the threads to finish, then accumulate the dBC_vec
    for(int t = 0; t < num_threads; ++t) {
        thread_vec[t].join();
    }
    for(int t = 0; t < num_threads; ++t) {
        for(node_id_t v = 0; v < dBC_vec.size(); ++v) {
            dBC_vec[v] += dBC_vec_vec[t][v];
        }
    }
    
    //now dBC_vec of each machine is ready.
    //master gets and accumulates all dBC_vec from everyone
    //printf("R[%d] done\n", rank);
//    if(rank == 0) {
//        //receive dBC_vec from everyone and accumulate
//        timer tm;
//        tm.start();
//        vector<double> rcv_dBC_vec;
//        rcv_dBC_vec.resize(dBC_vec.size());
//        for(int p = 1; p < size; ++p) {
//            MPI_Recv(&rcv_dBC_vec[0], rcv_dBC_vec.size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//            for(int i = 0; i < dBC_vec.size(); ++i) {
//                dBC_vec[i] += rcv_dBC_vec[i];
//            }
//        }
//        tm.stop();
//        printf("\t\tR[0] -- Accumulation time: [%f]\n", rank, tm.interval());
//    } else {
//        //just send dBC_vec
//        MPI_Send(&dBC_vec[0], dBC_vec.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//    }
    
    
    FILE* fout = fopen("dBC", "w");
    for(int i = 0; i < dBC_vec.size(); ++i) {
        fprintf(fout, "%f\n", dBC_vec[i]);
    }
    fclose(fout);
}

void iCentral_block(
        vector<double>*     dBC_vec,
        component_t*        comp,
        edge_t*             e,
        vector<node_id_t>*  source_vec,
        operation_t*        op
        )
{
    fill_vec<double>(*dBC_vec, comp->subgraph.size(), 0.0);
    iter_info_t iter_info;
    for(int i = 0; i < source_vec->size(); ++i) {
        node_id_t s = (*source_vec)[i];
        iCentral_iter(*dBC_vec, *comp, s, *e, iter_info, -1, false, *op);
    }
}


void iCentral_iter(
        vector<double>& dBC_vec,    // delta BC of vertices
        component_t&    comp,       // component could be BCC, MUC, or just a graph
        node_id_t       s,          // source of the iteration
        edge_t          e,          // inserted edge
        iter_info_t&    iter_info,  //TODO to be used later?
        int             dd,
        bool            use_d_1,
        operation_t      op
        )
{
    iter_info.init_all(comp.subgraph.size());
    if(dd > 0) {
        node_id_t tmp = e.first;
        e.first = e.second;
        e.second = tmp;
    }
    dd = abs(dd);
    if(op == INSERTION) {
        //NOTE: the d=1 optimization has a bug, however it was used in
        //      some of the performance results    
        //if(dd == 1 && use_d_1) {
        //    iter_info.init_new(comp.subgraph.size());//initializes the data structures required for BBFS_RBFS_d1
        //    BBFS_RBFS_d1(dBC_vec, iter_info, comp, s, e);
        //} else {
            BBFS(iter_info, comp, s);
            RBFS(dBC_vec, comp, s, iter_info, false, true);
            partial_BBFS_addition(iter_info, comp, s, e);
            RBFS(dBC_vec, comp, s, iter_info, true, false);
        //}
    } else if(op == DELETION) {
            BBFS(iter_info, comp, s);
            RBFS(dBC_vec, comp, s, iter_info, false, true);
            partial_BBFS_deletion(iter_info, comp, s, e);
            RBFS(dBC_vec, comp, s, iter_info, true, false);
    }

}


void brandes_iter(
        vector<double>& BC_vec,     // BC of vertices
        component_t&    comp,       // component could be BCC, MUC, or just a graph
        node_id_t       s,          // source of the iteration
        iter_info_t&    iter_info   //TODO to be used later?
        )
{
    iter_info.init_all(comp.subgraph.size());
    BBFS(iter_info, comp, s);
    RBFS(BC_vec, comp, s, iter_info, true, false);
}


void BBFS(
        iter_info_t&    iter_info,  // iteration info to be computed
        component_t&    comp,       // component
        node_id_t       s          // source of the iteration
        )
{
    subgraph_t& g = comp.subgraph;
    
    vector<vector<node_id_t> >&  P = iter_info.P; 
    vector<int>&                 sigma_vec = iter_info.sigma_vec;
    vector<int>&                 dist_vec = iter_info.dist_vec;
    vector<node_id_t>&           S = iter_info.S;
    
    queue<node_id_t>             Q;
    
    // Assumes iter_info is initialized
    sigma_vec[s] = 1;
    dist_vec[s] = 0;
    Q.push(s);
    while(!Q.empty()) {
        node_id_t v_i = Q.front(); Q.pop();
        S.push_back(v_i);
        for(int i = 0; i < g.nodes_vec[v_i].size(); ++i) {
            node_id_t v_n = g.nodes_vec[v_i][i];
            if(dist_vec[v_n] < 0) {
                Q.push(v_n);
                dist_vec[v_n] = dist_vec[v_i] + 1;
            }
            if(dist_vec[v_n] == dist_vec[v_i] + 1) {
                sigma_vec[v_n] = sigma_vec[v_n] + sigma_vec[v_i];
                P[v_n].push_back(v_i);
            }
        }
    }
}


void RBFS(
        vector<double>& dBC_vec,    // delta BC of vertices
        component_t&    comp,       // component could be BCC, MUC, or just a graph
        node_id_t       s,          // source of the iteration
        iter_info_t&    iter_info,  //
        bool            add,
        bool            sub
        )
{
    vector<vector<node_id_t> >&  P = iter_info.P; 
    vector<int>&                 sigma_vec = iter_info.sigma_vec;
    vector<double>&              delta_vec = iter_info.delta_vec;
    vector<double>&              Delta_vec = iter_info.Delta_vec;
    vector<node_id_t>&           S = iter_info.S;
    
    fill_vec<double>(delta_vec, comp.subgraph.size(), 0);
    
   
    if(comp.comp_type == GRAPH) {
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t w = S[i];
            for (int idx = 0; idx < P[w].size(); ++idx) {
                node_id_t v = P[w][idx];
                delta_vec[v] += (((double) sigma_vec[v] / sigma_vec[w]) * (1 + delta_vec[w]));
            }
            if (w != s) {
                if(add) dBC_vec[w] += delta_vec[w]/2.0;
                if(sub) dBC_vec[w] -= delta_vec[w]/2.0;
            }
        }
    } else if(comp.comp_type == BCC) { //BCC or MUC case: involved case, with external pairs contribution
        fill_vec<double>(Delta_vec, comp.subgraph.size(), 0);
        component_t::art_pt_map_t& art_pt_map = comp.art_pt_map;
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t v_n = S[i];
            if(   art_pt_map.find(s) != art_pt_map.end()
               && art_pt_map.find(v_n) != art_pt_map.end()
               ){
                int VG_s, VG_n;
                VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                VG_n = accumulate(art_pt_map[v_n].begin(), art_pt_map[v_n].end(), 0);
                int c_t = VG_s*VG_n;
                Delta_vec[v_n] =  Delta_vec[v_n] + c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)sigma_vec[v_p]/sigma_vec[v_n]);
                delta_vec[v_p] = delta_vec[v_p] + sp_sn*(1+delta_vec[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    Delta_vec[v_p] = Delta_vec[v_p] + Delta_vec[v_n]*sp_sn;
                }
            }
            if(s != v_n) {
                if(add) dBC_vec[v_n] += delta_vec[v_n]/2.0;
                if(sub) dBC_vec[v_n] -= delta_vec[v_n]/2.0;
            }
            if(art_pt_map.find(s) != art_pt_map.end()) {
                int VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                if(add) {
                    dBC_vec[v_n] += delta_vec[v_n]*VG_s;
                    dBC_vec[v_n] += Delta_vec[v_n]/2.0;
                }
                if(sub) {
                    dBC_vec[v_n] -= delta_vec[v_n]*VG_s;
                    dBC_vec[v_n] -= Delta_vec[v_n]/2.0;
                }
            }
        }
    }
    
}

void partial_BBFS_addition(
        iter_info_t&    iter_info,  // iteration info to be computed
        component_t&    comp,       // component
        node_id_t       s,          // source of the iteration
        edge_t          e           // inserted edge
        )
{
    subgraph_t& g = comp.subgraph;
    
    vector<vector<node_id_t> >&  P = iter_info.P; 
    vector<int>&                 path_cnt_vec = iter_info.sigma_vec;
    vector<int>&                 path_cnt_inc_map = iter_info.sigma_inc_vec;
    vector<int>&                 dist_map = iter_info.dist_vec;
    vector<bool>&                visited_vec = iter_info.visited_vec;
    vector<node_id_t>&           S = iter_info.S;
    
    queue<node_id_t>             Q;
    node_id_t src, dst;
    if(dist_map[e.first] > dist_map[e.second]) {
        src = e.second;
        dst = e.first;
    } else {
        src = e.first;
        dst = e.second;
    }
        
    //Compute new path counts and paths
    if(dist_map[dst] != dist_map[src]+1) {
        P[dst].clear();
        P[dst].push_back(src);
        path_cnt_vec[dst] = path_cnt_vec[src];
    } else {
        P[dst].push_back(src);
        path_cnt_vec[dst] += path_cnt_vec[src];
    }
        

    Q.push(dst);
    dist_map[dst] = dist_map[src] + 1;
    path_cnt_inc_map[dst] = path_cnt_vec[src];
    visited_vec[dst] = true;

    while (!Q.empty()) {
        node_id_t v = Q.front();
        Q.pop();
        vector<node_id_t> nbr_vec = g.get_nbrs(v);
        for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
            node_id_t w = nbr_vec[nbr_idx];
            if (dist_map[w] > (dist_map[v] + 1)) {
                dist_map[w] = dist_map[v] + 1;
                P[w].clear();
                P[w].push_back(v);
                path_cnt_vec[w] = 0;
                path_cnt_inc_map[w] = path_cnt_inc_map[v];
                path_cnt_vec[w] += path_cnt_inc_map[w];
                if (!visited_vec[w]) {
                    visited_vec[w] = true;
                    Q.push(w);
                }
            } else if (dist_map[w] == (dist_map[v] + 1)) {
                path_cnt_inc_map[w] += path_cnt_inc_map[v];
                path_cnt_vec[w] += path_cnt_inc_map[v];
                if (find(P[w].begin(), P[w].end(), v) == P[w].end()) {
                    P[w].push_back(v);
                }
                if (!visited_vec[w]) {
                    visited_vec[w] = true;
                    Q.push(w);
                }
            }
        }
    }


    //fix order of S
    //IMP::THIS CAN BE MADE MUCH BETTER!
    //HEAP FOR EXAMPLE
    //EVEN THE SWAPPING CAN BE DONE MORE EFFICIENTLY
    //for now it's not a bottleneck
    for (int i = 1; i < S.size(); ++i) {
        if (dist_map[S[i - 1]] > dist_map[S[i]]) {
            int j = i;
            while (dist_map[S[j - 1]] > dist_map[S[j]]) {
                node_id_t tmp = S[j - 1];
                S[j - 1] = S[j];
                S[j] = tmp;
                --j;
            }
        }
    }
}


/*
 TMP: this uses Brandes, full BFS, TODO: make efficient
 */
void partial_BBFS_deletion(
        iter_info_t&    iter_info,  // iteration info to be computed
        component_t&    comp,       // component
        node_id_t       s,          // source of the iteration
        edge_t          e           // deleted edge
        )
{
    subgraph_t& g = comp.subgraph;
    
    vector<vector<node_id_t> >&  P = iter_info.P; 
    vector<int>&                 sigma_vec = iter_info.sigma_vec;
    vector<int>&                 dist_vec = iter_info.dist_vec;
    vector<node_id_t>&           S = iter_info.S;
    
    queue<node_id_t>             Q;
    
    //for now everything is computed from scratch, so old values are irrelevant
    iter_info.init_all(g.size());
    // Assumes iter_info is initialized
    
    // IMP: careful bitch, graph_t is not thread safe!
    //g.remove_edge(e.first, e.second);
    sigma_vec[s] = 1;
    dist_vec[s] = 0;
    Q.push(s);
    while(!Q.empty()) {
        node_id_t v_i = Q.front(); Q.pop();
        S.push_back(v_i);
        for(int i = 0; i < g.nodes_vec[v_i].size(); ++i) {
            node_id_t v_n = g.nodes_vec[v_i][i];
            if(v_i == e.first && v_n == e.second || v_i == e.second && v_n == e.first) continue;
            if(dist_vec[v_n] < 0) {
                Q.push(v_n);
                dist_vec[v_n] = dist_vec[v_i] + 1;
            }
            if(dist_vec[v_n] == dist_vec[v_i] + 1) {
                sigma_vec[v_n] = sigma_vec[v_n] + sigma_vec[v_i];
                P[v_n].push_back(v_i);
            }
        }
    }
    //g.insert_edge(e.first, e.second);
}

void BBFS_RBFS_d1(
        vector<double>& dBC_vec,    // delta BC of vertices
        iter_info_t&    iter_info,  // iteration info to be computed
        component_t&    comp,       // component
        node_id_t       s,          // source of the iteration
        edge_t          e
        )
{
    subgraph_t& g = comp.subgraph;
    
    vector<vector<node_id_t> >&  P = iter_info.P; 
    vector<int>&                 sigma_vec = iter_info.sigma_vec;
    vector<int>&                 dist_vec = iter_info.dist_vec;
    vector<node_id_t>&           S = iter_info.S;
    vector<int>&                 new_sigma_vec = iter_info.new_sigma_vec;
    vector<double>&              delta_vec = iter_info.delta_vec;
    vector<double>&              Delta_vec = iter_info.Delta_vec;
    vector<double>&              new_delta_vec = iter_info.new_delta_vec;
    vector<double>&              new_Delta_vec = iter_info.new_Delta_vec;
    
    queue<node_id_t>             Q;
    
    // Assumes iter_info is initialized
    
    new_sigma_vec[s] = 1;//NEW
    sigma_vec[s] = 1;
    dist_vec[s] = 0;
    Q.push(s);
    while(!Q.empty()) {
        node_id_t v_i = Q.front(); Q.pop();
        S.push_back(v_i);
        if(v_i == e.second) {//NEW
            new_sigma_vec[v_i] += sigma_vec[e.first];
        }
        for(int i = 0; i < g.nodes_vec[v_i].size(); ++i) {
            node_id_t v_n = g.nodes_vec[v_i][i];
            if(dist_vec[v_n] < 0) {
                Q.push(v_n);
                dist_vec[v_n] = dist_vec[v_i] + 1;
            }
            if(dist_vec[v_n] == dist_vec[v_i] + 1) {
                sigma_vec[v_n] = sigma_vec[v_n] + sigma_vec[v_i];
                new_sigma_vec[v_n] += new_sigma_vec[v_i];//NEW
                P[v_n].push_back(v_i);
            }
        }
    }


    fill_vec<double>(delta_vec, comp.subgraph.size(), 0);
    fill_vec<double>(Delta_vec, comp.subgraph.size(), 0);
    
    fill_vec<double>(new_delta_vec, comp.subgraph.size(), 0);
    fill_vec<double>(new_Delta_vec, comp.subgraph.size(), 0);
    
    if(comp.comp_type == GRAPH) {
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t w = S[i];
            for (int idx = 0; idx < P[w].size(); ++idx) {
                node_id_t v = P[w][idx];
                delta_vec[v] += (((double) sigma_vec[v] / sigma_vec[w]) * (1 + delta_vec[w]));
                new_delta_vec[v] += (((double) new_sigma_vec[v] / new_sigma_vec[w]) * (1 + new_delta_vec[w]));
            }
            if(w == e.second) {
                node_id_t v_p = e.first;
                double new_sp_sn = ((double)new_delta_vec[v_p]/new_delta_vec[w]);
                new_delta_vec[v_p] = new_delta_vec[v_p] + new_sp_sn*(1+new_delta_vec[w]);
            }
            if (w != s) {
                dBC_vec[w] -= delta_vec[w]/2.0;
                dBC_vec[w] += new_delta_vec[w]/2.0;
            }
        }
    } else if(comp.comp_type == BCC) { //BCC or MUC case: involved case, with external pairs contribution
        component_t::art_pt_map_t& art_pt_map = comp.art_pt_map;
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t v_n = S[i];
            if(   art_pt_map.find(s) != art_pt_map.end()
               && art_pt_map.find(v_n) != art_pt_map.end()
               ){
                int VG_s, VG_n;
                VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                VG_n = accumulate(art_pt_map[v_n].begin(), art_pt_map[v_n].end(), 0);
                int c_t = VG_s*VG_n;
                Delta_vec[v_n] =  Delta_vec[v_n] + c_t;
                new_Delta_vec[v_n] =  new_Delta_vec[v_n] + c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)sigma_vec[v_p]/sigma_vec[v_n]);
                delta_vec[v_p] = delta_vec[v_p] + sp_sn*(1+delta_vec[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    Delta_vec[v_p] = Delta_vec[v_p] + Delta_vec[v_n]*sp_sn;
                }
                sp_sn = ((double)new_sigma_vec[v_p]/new_sigma_vec[v_n]);
                new_delta_vec[v_p] = new_delta_vec[v_p] + sp_sn*(1+new_delta_vec[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    new_Delta_vec[v_p] = new_Delta_vec[v_p] + new_Delta_vec[v_n]*sp_sn;
                }
            }
            //IMP: this is the only change that happens to P, @src should be added as parent for dst
            if(v_n == e.second) {
                node_id_t v_p = e.first;
                double new_sp_sn = ((double)new_delta_vec[v_p]/new_delta_vec[v_n]);
                new_delta_vec[v_p] = new_delta_vec[v_p] + new_sp_sn*(1+new_delta_vec[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    new_Delta_vec[v_p] = new_Delta_vec[v_p] + new_Delta_vec[v_n]*new_sp_sn;
                    dBC_vec[v_p] = dBC_vec[v_p] + new_Delta_vec[v_n]*new_sp_sn/2.0;
                }
            }
            if(s != v_n) {
                dBC_vec[v_n] -= delta_vec[v_n]/2.0;
                dBC_vec[v_n] += new_delta_vec[v_n]/2.0;
            }
            if(art_pt_map.find(s) != art_pt_map.end()) {
                int VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                dBC_vec[v_n] -= delta_vec[v_n]*VG_s;
                dBC_vec[v_n] -= Delta_vec[v_n]/2.0;
                dBC_vec[v_n] += new_delta_vec[v_n]*VG_s;
                dBC_vec[v_n] += new_Delta_vec[v_n]/2.0;
            }
        }
    }
}




/*
 * Parallel Brandes functions
 * Ziyad territory 
 */
void parallel_brandes(
        graph_t&        graph,
        vector<double>& BC_vec
        )
{
    //1. make a component (fill a subgraph)
    //2. loop through all sources
    BC_vec.resize(graph.size());
    fill(BC_vec.begin(), BC_vec.end(), 0.0);
    
    component_t comp;
    comp.comp_type = GRAPH;
    comp.subgraph.fill_graph(graph);
    
    // we want to do brandes_iter from each node in the graph
    // each thread will be responsible for doing brandes_iter from
    // a subset of the nodes in the graph
    //1. divide the nodes in the graph to vectors, each has a subset of nodes 
    //   a thread is responsible for
    //2. call brandes_block to get the contribution of the thread's nodes to the
    //   final BC value (fork and join)
    //3. accumulate the results in BC_vec

}

void brandes_block(
        vector<double>*     dBC_vec,
        component_t*        comp,
        vector<node_id_t>*  source_vec
        )
{
//    iter_info_t iter_info;
//    for(node_id_t s = 0; s < comp.subgraph.size(); ++s) {
//        brandes_iter(BC_vec, comp, s, iter_info);
//    }
}