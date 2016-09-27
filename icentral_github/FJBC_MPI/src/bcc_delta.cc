/* 
 * File:   graph_t.cc
 * Author: fuad
 * 
 * Created on January 6, 2014, 7:55 AM
 */
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include <stack>
#include <numeric>


#include "graph_t.h"
#include "utility.h"

using namespace std;

void bcc_delta_t::print()
{
    bcc_subgraph.print_graph(false);
    printf("Art. point and graph sizes:\n");
    for(art_pt_map_t::iterator it = art_pt_map.begin();
        it != art_pt_map.end();
        ++it) {
        printf("[%d]: ", it->first);
        for(int i = 0; i < it->second.size(); ++i) {
            printf("%d ", it->second[i]);
        }
        printf("\n");
    }
}

void graph_t::find_edge_bcc(bcc_delta_t&  bcc,
                                node_id_t src,
                                node_id_t dst)
{
    //a. find the bcc subgraph
    //b. find the articulation points in bcc
    //c. find the sizes of the subgraphs connected to each art point
    find_edge_bcc_subgraph(bcc.bcc_subgraph, src, dst);
    bcc.bcc_fast_subgraph.fill_graph(bcc.bcc_subgraph);
    vector<node_id_t> art_pt_vec;
    find_art_points(art_pt_vec);
    set<node_id_t> art_pt_set;
    for(int i = 0; i < art_pt_vec.size(); ++i) {
        art_pt_set.insert(art_pt_vec[i]);
    }
    
    graph_hash_t::nodes_map_t::iterator it;
    for(it = bcc.bcc_subgraph.nodes_map.begin();
        it != bcc.bcc_subgraph.nodes_map.end();
        ++it) {
        node_id_t v = it->first;
        if(art_pt_set.find(v) != art_pt_set.end()) {
            // the node is an articulation point and it belongs to the bcc
            // of interest, so it must be added to the art point of this bcc
            // along with the sizes of the subgraphs it connects the bcc to
            vector<int> subgraph_sz_vec;
            vector<bool> visited_vec;
            visited_vec.resize(size());
            fill(visited_vec.begin(), visited_vec.end(), false);
            visited_vec[v] = true;
            vector<node_id_t> v_nbr_vec = get_nbrs(v);
            for(int i = 0; i < v_nbr_vec.size(); ++i) {
                node_id_t u = v_nbr_vec[i];
                if(!bcc.bcc_subgraph.has_edge(v, u) && !visited_vec[u]) {
                    //do a dfs from this guy to figure out the size
                    //of the connected component connected to the bcc through v
                    //and this edge
                    int cnt = 0;
                    stack<node_id_t> S;
                    S.push(u);
                    visited_vec[u] = true;
                    while(!S.empty()) {
                        node_id_t vv = S.top();
                        S.pop();
                        cnt++;
                        vector<node_id_t> nbr_vec = get_nbrs(vv);
                        for(int i = 0; i < nbr_vec.size(); ++i) {
                            node_id_t nbr = nbr_vec[i];
                            if(!visited_vec[nbr]) {
                                visited_vec[nbr] = true;
                                S.push(nbr);
                            }
                        }
                    }
                    subgraph_sz_vec.push_back(cnt);
                }
            }
            bcc.art_pt_map.insert(make_pair(v, subgraph_sz_vec));
        }
    }
    
    //fix art_pt_map
    bcc_delta_t::art_pt_map_t new_art_pt_map;
    for(bcc_delta_t::art_pt_map_t::iterator it = bcc.art_pt_map.begin();
            it != bcc.art_pt_map.end();
            ++it) {
        node_id_t new_v = bcc.bcc_fast_subgraph.outin_label_map[it->first];
        new_art_pt_map.insert(make_pair(new_v, it->second));
    }
    bcc.art_pt_map = new_art_pt_map;
}

//////////////////////////////////////////////////
//IMP:: FOR NOW THIS GUY DOESN'T USE fast_subgraph
//////////////////////////////////////////////////
void bcc_delta_t::compute_bc(tr1_map_t(double)& bc_map,
                        node_id_t src,
                        node_id_t dst)
{

    //bc_map will have for each vertex in the muc it's betweenness centrality
    bcc_subgraph.remove_edge(src, dst);
    tr1_map_t(int) d_src_vec, d_dst_vec;
    bcc_subgraph.find_sssp(src, d_src_vec);
    bcc_subgraph.find_sssp(dst, d_dst_vec);
    
    //this must be commented in general.. if it's there it makes bc_map contain deltas
    //bcc_subgraph.i_fill_map<double>(bc_map, 0);
    
    for(graph_hash_t::nodes_map_t::iterator
                        it =  bcc_subgraph.nodes_map.begin();
                        it != bcc_subgraph.nodes_map.end();
                        ++it) {
        node_id_t s = it->first;
        if(d_src_vec[s] != d_dst_vec[s]) {
            i_iteration(s, 
                    src, 
                    dst, 
                    d_src_vec[s], 
                    d_dst_vec[s], 
                    bc_map);
        }
    }
    bcc_subgraph.insert_edge(src, dst);
}

void bcc_delta_t::compute_bc_exp(tr1_map_t(double)& bc_map,
                        node_id_t src,
                        node_id_t dst,
                        bcc_stat_t& bcc_stat,
                        int max_iter_d1,
                        int max_iter_d2)
{
    src = bcc_fast_subgraph.outin_label_map[src];
    dst = bcc_fast_subgraph.outin_label_map[dst];
    
    timer tm;
    //bc_map will have for each vertex in the muc it's betweenness centrality
    bcc_fast_subgraph.remove_edge(src, dst);
    vector<int> d_src_vec, d_dst_vec;
    tm.start();
    bcc_fast_subgraph.find_sssp(src, d_src_vec);
    bcc_fast_subgraph.find_sssp(dst, d_dst_vec);
    tm.stop();
    bcc_stat.sssp_tm = tm.interval();
    
    //this must be commented in general.. if it's there it makes bc_map contain deltas
    //bcc_subgraph.i_fill_map<double>(bc_map, 0);
    
    int    cnt_arr[] = {0,0,0};
    double tot_arr[] = {0,0,0};
    for(int i = 0; i < bcc_fast_subgraph.nodes_vec.size(); ++i) {
        node_id_t s = i;
        int diff = d_src_vec[s] - d_dst_vec[s];
        int abs_diff = abs(diff);
        tm.start();
        if(d_src_vec[s] != d_dst_vec[s]) {
            i_iteration(s, src, dst, d_src_vec[s], d_dst_vec[s], bc_map);
        }
        tm.stop();
        switch(abs_diff) {
            case 0:
                cnt_arr[0]++;
                tot_arr[0]+=tm.interval();
                break;
            case 1:
                cnt_arr[1]++;
                tot_arr[1]+=tm.interval();
                break;
            default:
                cnt_arr[2]++;
                tot_arr[2]+=tm.interval();
        }
        if(tot_arr[1] == max_iter_d1 && tot_arr[2] == max_iter_d2)
            break;
    }
    bcc_fast_subgraph.insert_edge(src, dst);
    
    bcc_stat.num_d0_iter = cnt_arr[0];
    bcc_stat.num_d1_iter = cnt_arr[1];
    bcc_stat.num_d2_iter = cnt_arr[2];
    bcc_stat.tot_d0_tm   = tot_arr[0];
    bcc_stat.tot_d1_tm   = tot_arr[1];
    bcc_stat.tot_d2_tm   = tot_arr[2];
}

void bcc_delta_t::compute_bc_maxiter_exp(tr1_map_t(double)& bc_map,
                        node_id_t src,
                        node_id_t dst,
                        bcc_stat_t& bcc_stat,
                        int max_iter_d1,
                        int max_iter_d2)
{
    src = bcc_fast_subgraph.outin_label_map[src];
    dst = bcc_fast_subgraph.outin_label_map[dst];
    
    timer tm;
    //bc_map will have for each vertex in the muc it's betweenness centrality
    bcc_fast_subgraph.remove_edge(src, dst);
    vector<int> d_src_vec, d_dst_vec;
    tm.start();
    bcc_fast_subgraph.find_sssp(src, d_src_vec);
    bcc_fast_subgraph.find_sssp(dst, d_dst_vec);
    tm.stop();
    bcc_stat.sssp_tm = tm.interval();
    
    //this must be commented in general.. if it's there it makes bc_map contain deltas
    //bcc_subgraph.i_fill_map<double>(bc_map, 0);
    vector<node_id_t> d1_s_vec, d2_s_vec;
    
    int    cnt_arr[] = {0,0,0};
    double tot_arr[] = {0,0,0};
    for(int i = 0; i < bcc_fast_subgraph.nodes_vec.size(); ++i) {
        node_id_t s = i;
        int diff = d_src_vec[s] - d_dst_vec[s];
        int abs_diff = abs(diff);        
        switch(abs_diff) {
            case 0:  cnt_arr[0]++; break;
            case 1:  d1_s_vec.push_back(s); break;
            default: d2_s_vec.push_back(s);
        }
    }
    
    int d1_num_iter = min(max_iter_d1, (int)d1_s_vec.size());
    int d2_num_iter = min(max_iter_d2, (int)d2_s_vec.size());
    for(int i = 0; i < d1_num_iter; ++i) {
        node_id_t s = d1_s_vec[i];
        tm.start();
        i_iteration(s, src, dst, d_src_vec[s], d_dst_vec[s], bc_map);
        tm.stop();
        cnt_arr[1]++;
        tot_arr[1]+=tm.interval();
    }
    for(int i = 0; i < d2_num_iter; ++i) {
        node_id_t s = d2_s_vec[i];
        tm.start();
        i_iteration(s, src, dst, d_src_vec[s], d_dst_vec[s], bc_map);
        tm.stop();
        cnt_arr[2]++;
        tot_arr[2]+=tm.interval();
    }
    bcc_fast_subgraph.insert_edge(src, dst);
    
    bcc_stat.num_d0_iter = cnt_arr[0];
    bcc_stat.num_d1_iter = cnt_arr[1];
    bcc_stat.num_d2_iter = cnt_arr[2];
    bcc_stat.tot_d0_tm   = tot_arr[0];
    bcc_stat.tot_d1_tm   = tot_arr[1];
    bcc_stat.tot_d2_tm   = tot_arr[2];
}

void bcc_delta_t::i_iteration(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_vec)
{
    //make sure that @src is the closer to source node
    if(d_src > d_dst) {
        node_id_t       tmp_v;
        int             tmp_i;
        tmp_v = src; src = dst; dst = tmp_v;
        tmp_i = d_src; d_src = d_dst; d_dst = tmp_i;
    }
    
    if(d_dst - d_src == 1) {
        //dbg_iteration(s, src, dst, d_src, d_dst, bc_vec);
        i_iteration_1(s, src, dst, d_src, d_dst, bc_vec);
        //i_iteration_2(s, src, dst, d_src, d_dst, bc_vec);
    } else {
        //dbg_iteration(s, src, dst, d_src, d_dst, bc_vec);
        i_iteration_2(s, src, dst, d_src, d_dst, bc_vec);
    }
}

void bcc_delta_t::dbg_iteration(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_map)
{
        tr1_map_t(vector<node_id_t> )&  P = bcc_subgraph.P; 
        tr1_map_t(int)&                 path_cnt_map = bcc_subgraph.path_cnt_vec;
        tr1_map_t(int)&                 dist_map = bcc_subgraph.dist_vec;
        tr1_map_t(double)&              pair_dep_map = bcc_subgraph.pair_dep_vec;
        tr1_map_t(double)&              sigma_t_map = bcc_subgraph.sigma_t_map;
        queue<node_id_t>&               Q = bcc_subgraph.Q;
        vector<node_id_t>               S;
        
        bcc_subgraph.init_maps();           
        path_cnt_map[s] = 1;
        dist_map[s] = 0;
        Q.push(s);
        while(!Q.empty()) {
            node_id_t v_i = Q.front(); Q.pop();
            S.push_back(v_i);
            for(int i = 0; i < bcc_subgraph.nodes_map[v_i].nbrs_vec.size(); ++i) {
                node_id_t v_n = bcc_subgraph.nodes_map[v_i].nbrs_vec[i];
                if(dist_map[v_n] < 0) {
                    Q.push(v_n);
                    dist_map[v_n] = dist_map[v_i] + 1;
                }
                if(dist_map[v_n] == dist_map[v_i] + 1) {
                    path_cnt_map[v_n] = path_cnt_map[v_n] + path_cnt_map[v_i];
                    P[v_n].push_back(v_i);
                }
            }
        }
        //RBFS to subtract old pair dependency
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t v_n = S[i];
            if(   art_pt_map.find(s) != art_pt_map.end()
               && art_pt_map.find(v_n) != art_pt_map.end()
               && s != v_n) {
                int VG_s, VG_n;
                VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                VG_n = accumulate(art_pt_map[v_n].begin(), art_pt_map[v_n].end(), 0);
                int c_t = VG_s*VG_n;
                sigma_t_map[v_n] =  sigma_t_map[v_n] + c_t;
                /*this guy must not change!*/
                //bc_map[v_n] = bc_map[v_n] - c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)path_cnt_map[v_p]/path_cnt_map[v_n]);
                pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn*(1+pair_dep_map[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
                    bc_map[v_p] = bc_map[v_p] - sigma_t_map[v_n]*sp_sn/2.0;
                }
            }
            if(s != v_n) {
                bc_map[v_n] = bc_map[v_n] - pair_dep_map[v_n]/2.0;
            }
            if(art_pt_map.find(s) != art_pt_map.end()) {
                int VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                bc_map[v_n] = bc_map[v_n] - pair_dep_map[v_n]*VG_s;
            }
        }
        
        bcc_subgraph.insert_edge(src, dst);
        bcc_subgraph.init_maps();
        path_cnt_map[s] = 1;
        dist_map[s] = 0;
        Q.push(s);
        S.clear();
        while(!Q.empty()) {
            node_id_t v_i = Q.front(); Q.pop();
            S.push_back(v_i);
            for(int i = 0; i < bcc_subgraph.nodes_map[v_i].nbrs_vec.size(); ++i) {
                node_id_t v_n = bcc_subgraph.nodes_map[v_i].nbrs_vec[i];
                if(dist_map[v_n] < 0) {
                    Q.push(v_n);
                    dist_map[v_n] = dist_map[v_i] + 1;
                }
                if(dist_map[v_n] == dist_map[v_i] + 1) {
                    path_cnt_map[v_n] = path_cnt_map[v_n] + path_cnt_map[v_i];
                    P[v_n].push_back(v_i);
                }
            }
        }
        
        //DBG
//        printf("===========================\nSource: %d\n", s);
//        for(tr1_map_t(double)::iterator it = pair_dep_map.begin();
//                it != pair_dep_map.end();
//                ++it) {
//            printf("pd[%d]: %f           ", it->first, it->second);
//            printf("# parents: %d\n", P[it->first].size());
//        }
        //////
        
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t v_n = S[i];
            if(   art_pt_map.find(s) != art_pt_map.end()
               && art_pt_map.find(v_n) != art_pt_map.end()
               && s != v_n) {
                int VG_s, VG_n;
                VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                VG_n = accumulate(art_pt_map[v_n].begin(), art_pt_map[v_n].end(), 0);
                int c_t = VG_s*VG_n;
                sigma_t_map[v_n] =  sigma_t_map[v_n] + c_t;
                /*this guy must not change!*/
                //bc_map[v_n] = bc_map[v_n] + c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)path_cnt_map[v_p]/path_cnt_map[v_n]);
                pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn*(1+pair_dep_map[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
                    bc_map[v_p] = bc_map[v_p] + sigma_t_map[v_n]*sp_sn/2.0;
                }
            }
            if(s != v_n) {
                bc_map[v_n] = bc_map[v_n] + pair_dep_map[v_n]/2.0;
            }
            if(art_pt_map.find(s) != art_pt_map.end()) {
                int VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                bc_map[v_n] = bc_map[v_n] + pair_dep_map[v_n]*VG_s;
            }
        }
        
        //DBG
//        printf("===========================\nSource: %d\n", s);
//        for(tr1_map_t(int)::iterator it = path_cnt_map.begin();
//                it != path_cnt_map.end();
//                ++it) {
//            printf("cnt[%d]: %d           ", it->first, it->second);
//            printf("# parents: %d\n", P[it->first].size());
//        }
        /////
        
        //DBG
//        printf("========================\nSource: %d\n", s);
//        for(int i = S.size()-1; i >= 0; --i) {
//            printf("%d\n", S[i]);
//        }
        //////
        
        //DBG print all:
//        printf("===========================\nSource: %d\n", s);
//        for(tr1_map_t(int)::iterator it = path_cnt_map.begin();
//                it != path_cnt_map.end();
//                ++it) {
//            printf("cnt[%d]: %d\n", it->first, it->second);
//            for(int i = 0; i < P[it->first].size(); ++i) {
//                printf("        P[%d][%d]: %d\n", it->first, i, P[it->first][i]);
//            }
//        }
//        printf("S:  ");
//        for(int i = S.size()-1; i >= 0; --i) {
//            printf("%d  ", S[i]);
//        }
//        printf("\n");
        ///////////
        

        bcc_subgraph.remove_edge(src, dst);
}

template <typename T>
void fill_vec(vector<T>& vec, T val, int size)
{
    vec.resize(size);
    fill(vec.begin(), vec.end(), val);
}

void bcc_delta_t::i_iteration_1(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_map)
{
        vector<vector<node_id_t> >&  P = bcc_fast_subgraph.P; 
        vector<int>&                 path_cnt_map = bcc_fast_subgraph.path_cnt_vec;
        vector<int>&                 new_path_cnt_map = bcc_fast_subgraph.new_path_cnt_vec;
        vector<int>&                 dist_map = bcc_fast_subgraph.dist_vec;
        vector<double>&              pair_dep_map = bcc_fast_subgraph.pair_dep_vec;
        vector<double>&              new_pair_dep_map = bcc_fast_subgraph.new_pair_dep_vec;
        vector<double>&              sigma_t_map = bcc_fast_subgraph.sigma_t_map;
        vector<double>&              new_sigma_t_map = bcc_fast_subgraph.new_sigma_t_map; 
        queue<node_id_t>&            Q = bcc_fast_subgraph.Q;
        stack<node_id_t>             S;
        bcc_fast_subgraph.init_maps();
        
//    int sz = bcc_fast_subgraph.nodes_vec.size();
//        vector<vector<node_id_t> >  P(sz, vector<node_id_t>());// = bcc_fast_subgraph.P; 
//        vector<int>                 path_cnt_map(sz, 0);// = bcc_fast_subgraph.path_cnt_vec;
//        vector<int>                 new_path_cnt_map(sz, 0);// = bcc_fast_subgraph.new_path_cnt_vec;
//        vector<int>                 dist_map(sz, -1);// = bcc_fast_subgraph.dist_vec;
//        vector<double>              pair_dep_map(sz, 0);// = bcc_fast_subgraph.pair_dep_vec;
//        vector<double>              new_pair_dep_map(sz, 0);// = bcc_fast_subgraph.new_pair_dep_vec;
//        vector<double>              sigma_t_map(sz, 0);// = bcc_fast_subgraph.sigma_t_map;
//        vector<double>              new_sigma_t_map(sz, 0);// = bcc_fast_subgraph.new_sigma_t_map; 
//        queue<node_id_t>            Q;// = bcc_fast_subgraph.Q;
//        stack<node_id_t>            S;
        
        
        path_cnt_map[s] = 1;
        new_path_cnt_map[s] = 1;//NEW
        dist_map[s] = 0;
        
        Q.push(s);
        while(!Q.empty()) {
            node_id_t v_i = Q.front(); Q.pop();
            S.push(v_i);
            if(v_i == dst) {//NEW
                new_path_cnt_map[v_i] += path_cnt_map[src];
            }
            for(int i = 0; i < bcc_fast_subgraph.nodes_vec[v_i].size(); ++i) {
                node_id_t v_n = bcc_fast_subgraph.nodes_vec[v_i][i];
                if(dist_map[v_n] < 0) {
                    Q.push(v_n);
                    dist_map[v_n] = dist_map[v_i] + 1;
                }
                if(dist_map[v_n] == dist_map[v_i] + 1) {
                    path_cnt_map[v_n] = path_cnt_map[v_n] + path_cnt_map[v_i];
                    new_path_cnt_map[v_n] += new_path_cnt_map[v_i];//NEW
                    P[v_n].push_back(v_i);
                }
            }
        }
        
        while(!S.empty()) {
            node_id_t v_n = S.top(); S.pop();
            if(   art_pt_map.find(s) != art_pt_map.end()
               && art_pt_map.find(v_n) != art_pt_map.end()
               && s != v_n) {
                int VG_s, VG_n;
                VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                VG_n = accumulate(art_pt_map[v_n].begin(), art_pt_map[v_n].end(), 0);
                int c_t = VG_s*VG_n;
                sigma_t_map[v_n] =  sigma_t_map[v_n] + c_t;
                new_sigma_t_map[v_n] = new_sigma_t_map[v_n] + c_t;//NEW
                /*this guy must not change!*/
                //bc_map[v_n] = bc_map[v_n] - c_t;
                //bc_map[v_n] = bc_map[v_n] + c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)path_cnt_map[v_p]/path_cnt_map[v_n]);
                pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn*(1+pair_dep_map[v_n]);
                double new_sp_sn = ((double)new_path_cnt_map[v_p]/new_path_cnt_map[v_n]);
                new_pair_dep_map[v_p] = new_pair_dep_map[v_p] + new_sp_sn*(1+new_pair_dep_map[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
                    new_sigma_t_map[v_p] = new_sigma_t_map[v_p] + new_sigma_t_map[v_n]*new_sp_sn;
                    
                    node_id_t new_v_p = bcc_fast_subgraph.inout_label_map[v_p];
                    bc_map[new_v_p] = bc_map[new_v_p] - sigma_t_map[v_n]*sp_sn/2.0;
                    bc_map[new_v_p] = bc_map[new_v_p] + new_sigma_t_map[v_n]*new_sp_sn/2.0;
                }
            }
            //IMP: this is the only change that happens to P, @src should be added as parent for dst
            if(v_n == dst) {
                node_id_t v_p = src;
                double new_sp_sn = ((double)new_path_cnt_map[v_p]/new_path_cnt_map[v_n]);
                new_pair_dep_map[v_p] = new_pair_dep_map[v_p] + new_sp_sn*(1+new_pair_dep_map[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    new_sigma_t_map[v_p] = new_sigma_t_map[v_p] + new_sigma_t_map[v_n]*new_sp_sn;
                    
                    node_id_t new_v_p = bcc_fast_subgraph.inout_label_map[v_p];
                    bc_map[new_v_p] = bc_map[new_v_p] + new_sigma_t_map[v_n]*new_sp_sn/2.0;
                }
            }
            if(s != v_n) {
                node_id_t new_v_n = bcc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] - pair_dep_map[v_n]/2.0;
                bc_map[new_v_n] = bc_map[new_v_n] + new_pair_dep_map[v_n]/2.0;
            }
            if(art_pt_map.find(s) != art_pt_map.end()) {
                int VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                
                node_id_t new_v_n = bcc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] - pair_dep_map[v_n]*VG_s;
                bc_map[new_v_n] = bc_map[new_v_n] + new_pair_dep_map[v_n]*VG_s;
            }
        }
}

void bcc_delta_t::i_iteration_2(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_map)
{
        vector<vector<node_id_t> >&  P = bcc_fast_subgraph.P; 
        vector<int>&                 path_cnt_map = bcc_fast_subgraph.path_cnt_vec;
        vector<int>&                 path_cnt_inc_map = bcc_fast_subgraph.path_cnt_inc_vec;
        vector<int>&                 dist_map = bcc_fast_subgraph.dist_vec;
        vector<double>&              pair_dep_map = bcc_fast_subgraph.pair_dep_vec;
        vector<double>&              sigma_t_map = bcc_fast_subgraph.sigma_t_map;
        vector<bool>&                visited_vec = bcc_fast_subgraph.visited_vec;
        queue<node_id_t>&            Q = bcc_fast_subgraph.Q;
        vector<node_id_t>            S;
        
        bcc_fast_subgraph.init_maps();           

        path_cnt_map[s] = 1;
        dist_map[s] = 0;
        Q.push(s);
        while(!Q.empty()) {
            node_id_t v_i = Q.front(); Q.pop();
            S.push_back(v_i);
            for(int i = 0; i < bcc_fast_subgraph.nodes_vec[v_i].size(); ++i) {
                node_id_t v_n = bcc_fast_subgraph.nodes_vec[v_i][i];
                if(dist_map[v_n] < 0) {
                    Q.push(v_n);
                    dist_map[v_n] = dist_map[v_i] + 1;
                }
                if(dist_map[v_n] == dist_map[v_i] + 1) {
                    path_cnt_map[v_n] = path_cnt_map[v_n] + path_cnt_map[v_i];
                    P[v_n].push_back(v_i);
                }
            }
        }
        
        /*
         * steps:
         * 1. do the reverse BFS and subtract the olds
         * 2. compute the new counts
         * 3. fix the order of S
         * 4. add the new increments
         */
        //RBFS to subtract old pair dependency
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t v_n = S[i];
            if(   art_pt_map.find(s) != art_pt_map.end()
               && art_pt_map.find(v_n) != art_pt_map.end()
               ){//&& s != v_n) {
                int VG_s, VG_n;
                VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                VG_n = accumulate(art_pt_map[v_n].begin(), art_pt_map[v_n].end(), 0);
                int c_t = VG_s*VG_n;
                sigma_t_map[v_n] =  sigma_t_map[v_n] + c_t;
                /*this guy must not change!*/
                //bc_map[v_n] = bc_map[v_n] - c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)path_cnt_map[v_p]/path_cnt_map[v_n]);
                pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn*(1+pair_dep_map[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
                    
//                    node_id_t new_v_p = bcc_fast_subgraph.inout_label_map[v_p];
//                    bc_map[new_v_p] = bc_map[new_v_p] - sigma_t_map[v_n]*sp_sn/2.0;
                }
            }
            if(s != v_n) {
                node_id_t new_v_n = bcc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] - pair_dep_map[v_n]/2.0;
            }
            if(art_pt_map.find(s) != art_pt_map.end()) {
                int VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                
                node_id_t new_v_n = bcc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] - pair_dep_map[v_n]*VG_s;
                
                bc_map[new_v_n] = bc_map[new_v_n] - sigma_t_map[v_n]/2.0;
            }
        }
                
        //Compute new path counts and paths
        if(dist_map[dst] != dist_map[src]+1) {
            P[dst].clear();
            P[dst].push_back(src);
            path_cnt_map[dst] = path_cnt_map[src];
        } else {
            P[dst].push_back(src);
            path_cnt_map[dst] += path_cnt_map[src];
        }
        
        
        Q.push(dst);
        dist_map[dst] = dist_map[src] + 1;
        //P[dst].clear();
        //P[dst].push_back(src);
        path_cnt_inc_map[dst] = path_cnt_map[src];
        //path_cnt_map[dst] = path_cnt_inc_map[dst];
        visited_vec[dst] = true;
        
        while (!Q.empty()) {
            node_id_t v = Q.front();
            Q.pop();
            vector<node_id_t> nbr_vec = bcc_fast_subgraph.get_nbrs(v);
            for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
                node_id_t w = nbr_vec[nbr_idx];
                if (dist_map[w] > (dist_map[v] + 1)) {
                    dist_map[w] = dist_map[v] + 1;
                    P[w].clear();
                    P[w].push_back(v);
                    path_cnt_map[w] = 0;
                    path_cnt_inc_map[w] = path_cnt_inc_map[v];
                    path_cnt_map[w] += path_cnt_inc_map[w];
                    if(!visited_vec[w]) {
                        visited_vec[w] = true;
                        Q.push(w);
                    }
                } else if (dist_map[w] == (dist_map[v] + 1)) {
                    path_cnt_inc_map[w] += path_cnt_inc_map[v];
                    path_cnt_map[w] += path_cnt_inc_map[v];
                    if(find(P[w].begin(), P[w].end(), v) == P[w].end()) {
                        P[w].push_back(v);
                    }
                    if(!visited_vec[w]) {
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
        for(int i = 1; i < S.size(); ++i) {
            if(dist_map[S[i-1]] > dist_map[S[i]]) {
                int j = i;
                while(dist_map[S[j-1]] > dist_map[S[j]]) {
                    node_id_t tmp = S[j-1];
                    S[j-1] = S[j];
                    S[j] = tmp;
                    --j;
                }
            }
        }

        //DBG
//        printf("===========================\nSource: %d\n", s);
//        for(tr1_map_t(int)::iterator it = path_cnt_map.begin();
//                it != path_cnt_map.end();
//                ++it) {
//            printf("cnt[%d]: %d           ", it->first, it->second);
//            printf("# parents: %d\n", P[it->first].size());
//        }
        /////        
        //DBG
//        printf("========================\nSource: %d\n", s);
//        for(int i = S.size()-1; i >= 0; --i) {
//            printf("%d\n", S[i]);
//        }
        //////
        

        //DBG print all:
//        printf("===========================\nSource: %d\n", s);
//        for(tr1_map_t(int)::iterator it = path_cnt_map.begin();
//                it != path_cnt_map.end();
//                ++it) {
//            printf("cnt[%d]: %d\n", it->first, it->second);
//            for(int i = 0; i < P[it->first].size(); ++i) {
//                printf("        P[%d][%d]: %d\n", it->first, i, P[it->first][i]);
//            }
//        }
//        printf("S:  ");
//        for(int i = S.size()-1; i >= 0; --i) {
//            printf("%d  ", S[i]);
//        }
//        printf("\n");
        ///////////
        

        //DBG
//        printf("===========================\nSource: %d\n", s);
//        for(tr1_map_t(double)::iterator it = pair_dep_map.begin();
//                it != pair_dep_map.end();
//                ++it) {
//            printf("pd[%d]: %f           ", it->first, it->second);
//            printf("# parents: %d\n", P[it->first].size());
//        }
        //////
        
        bcc_fast_subgraph.i_fill_map<double>(pair_dep_map, 0);
        bcc_fast_subgraph.i_fill_map<double>(sigma_t_map, 0);
        //RBFS to add the new pair dependencies
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t v_n = S[i];
            if(   art_pt_map.find(s) != art_pt_map.end()
               && art_pt_map.find(v_n) != art_pt_map.end()
               ){//&& s != v_n) {
                int VG_s, VG_n;
                VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                VG_n = accumulate(art_pt_map[v_n].begin(), art_pt_map[v_n].end(), 0);
                int c_t = VG_s*VG_n;
                sigma_t_map[v_n] =  sigma_t_map[v_n] + c_t;
                /*this guy must not change!*/
                //bc_map[v_n] = bc_map[v_n] + c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)path_cnt_map[v_p]/path_cnt_map[v_n]);
                pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn*(1+pair_dep_map[v_n]);
                if(art_pt_map.find(s) != art_pt_map.end()) {
                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
                    
//                    node_id_t new_v_p = bcc_fast_subgraph.inout_label_map[v_p];
//                    bc_map[new_v_p] = bc_map[new_v_p] + sigma_t_map[v_n]*sp_sn/2.0;
                }
            }
            if(s != v_n) {
                node_id_t new_v_n = bcc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] + pair_dep_map[v_n]/2.0;
            }
            if(art_pt_map.find(s) != art_pt_map.end()) {
                int VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
                
                node_id_t new_v_n = bcc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] + pair_dep_map[v_n]*VG_s;
                
                bc_map[new_v_n] = bc_map[new_v_n] + sigma_t_map[v_n]/2.0;
            }
        }

}