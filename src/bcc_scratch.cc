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

using namespace std;

void bcc_scratch_t::print()
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

void graph_t::find_edge_bcc(bcc_scratch_t& bcc,
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



void bcc_scratch_t::compute_bc(tr1::unordered_map<node_id_t, double>& bc_map, int max_iter)
{
    //bc_map will have for each vertex in the muc it's betweenness centrality
    //init bc of all nodes to zero
    //XXX why fill this bc_map? Should have all vertices in the graph, here we increment/decrement
    bcc_subgraph.i_fill_map<double>(bc_map, 0);
    
    //do BFS's from the nodes in the bcc
    if(max_iter == -1) {
        for(int i = 0; i < bcc_fast_subgraph.nodes_vec.size(); ++i) {
            node_id_t s = i;
            bc_iter(s, bc_map);
        }  
    } else {
        int cnt = 0;
        for(int i = 0; i < bcc_fast_subgraph.nodes_vec.size(); ++i) {
            node_id_t s = i;
            ++cnt;
            bc_iter(s, bc_map);
            if(cnt == max_iter)
                break;
        }  
    }

    if(max_iter == -1) {
        for(art_pt_map_t::iterator
            it =  art_pt_map.begin();
            it != art_pt_map.end();
            ++it) {
            vector<int> size_vec;
            size_vec = it->second;
            if(size_vec.size() > 1) {
                int sub = 0;
                for(int i = 0; i < size_vec.size(); ++i) {
                    sub += size_vec[i]*size_vec[i];
                }
                int VG_i = accumulate(it->second.begin(), it->second.end(), 0);
                bc_map[it->first] = bc_map[it->first] + VG_i*VG_i - sub;
            }
        }
        for (tr1::unordered_map<node_id_t, double>::iterator
            it = bc_map.begin();
                it != bc_map.end();
                ++it) {
            it->second = it->second / 2;
        }
    }
}


void bcc_scratch_t::bc_iter(node_id_t s, tr1::unordered_map<node_id_t,double>& bc_map) {
    vector<vector<node_id_t> >&  P = bcc_fast_subgraph.P;
    vector<int>&                 path_cnt_map = bcc_fast_subgraph.path_cnt_vec;
    vector<int>&                 dist_map = bcc_fast_subgraph.dist_vec;
    vector<double>&              pair_dep_map = bcc_fast_subgraph.pair_dep_vec;
    vector<double>&              sigma_t_map = bcc_fast_subgraph.sigma_t_map;
    queue<node_id_t>&            Q = bcc_fast_subgraph.Q;
    stack<node_id_t>             S;

    bcc_fast_subgraph.init_maps();
    path_cnt_map[s] = 1;
    dist_map[s] = 0;

    Q.push(s);
    while (!Q.empty()) {
        node_id_t v_i = Q.front();
        Q.pop();
        S.push(v_i);
        for(int i = 0; i < bcc_fast_subgraph.nodes_vec[v_i].size(); ++i) {
            node_id_t v_n = bcc_fast_subgraph.nodes_vec[v_i][i];
            if (dist_map[v_n] < 0) {
                Q.push(v_n);
                dist_map[v_n] = dist_map[v_i] + 1;
            }
            if (dist_map[v_n] == dist_map[v_i] + 1) {
                path_cnt_map[v_n] = path_cnt_map[v_n] + path_cnt_map[v_i];
                P[v_n].push_back(v_i);
            }
        }
    }
    while (!S.empty()) {
        node_id_t v_n = S.top();
        S.pop();
        if (art_pt_map.find(s) != art_pt_map.end()
                && art_pt_map.find(v_n) != art_pt_map.end()
                && s != v_n) {
            int VG_s, VG_n;
            VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
            VG_n = accumulate(art_pt_map[v_n].begin(), art_pt_map[v_n].end(), 0);
            int c_t = VG_s*VG_n;
            sigma_t_map[v_n] = sigma_t_map[v_n] + c_t;
            
            node_id_t new_v_n = bcc_fast_subgraph.inout_label_map[v_n];
            bc_map[new_v_n] = bc_map[new_v_n] + c_t;
        }
        for (int i = 0; i < P[v_n].size(); ++i) {
            node_id_t v_p = P[v_n][i];
            double sp_sn = ((double) path_cnt_map[v_p] / path_cnt_map[v_n]);
            pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn * (1 + pair_dep_map[v_n]);
            if (art_pt_map.find(s) != art_pt_map.end()) {
                sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n] * sp_sn;
                
                node_id_t new_v_p = bcc_fast_subgraph.inout_label_map[v_p];
                bc_map[new_v_p] = bc_map[new_v_p] + sigma_t_map[v_n] * sp_sn;
            }
        }
        if (s != v_n) {
            node_id_t new_v_n = bcc_fast_subgraph.inout_label_map[v_n];
            bc_map[new_v_n] = bc_map[new_v_n] + pair_dep_map[v_n];
        }
        if (art_pt_map.find(s) != art_pt_map.end()) {
            int VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
            
            node_id_t new_v_n = bcc_fast_subgraph.inout_label_map[v_n];
            bc_map[new_v_n] = bc_map[new_v_n] + pair_dep_map[v_n] * VG_s * 2;
        }
    }    
}

///* 
// * File:   graph_t.cc
// * Author: fuad
// * 
// * Created on January 6, 2014, 7:55 AM
// */
//#include <iostream>
//#include <fstream>
//#include <stdio.h>
//#include <algorithm>
//#include <queue>
//#include <stack>
//#include <numeric>
//
//
//#include "graph_t.h"
//
//using namespace std;
//
//void bcc_scratch_t::print()
//{
//    bcc_subgraph.print_graph(false);
//    printf("Art. point and graph sizes:\n");
//    for(art_pt_map_t::iterator it = art_pt_map.begin();
//        it != art_pt_map.end();
//        ++it) {
//        printf("[%d]: ", it->first);
//        for(int i = 0; i < it->second.size(); ++i) {
//            printf("%d ", it->second[i]);
//        }
//        printf("\n");
//    }
//}
//
//void graph_t::find_edge_bcc(bcc_scratch_t& bcc,
//                                node_id_t src,
//                                node_id_t dst)
//{
//    //a. find the bcc subgraph
//    //b. find the articulation points in bcc
//    //c. find the sizes of the subgraphs connected to each art point
//    find_edge_bcc_subgraph(bcc.bcc_subgraph, src, dst);
//    vector<node_id_t> art_pt_vec;
//    find_art_points(art_pt_vec);
//    set<node_id_t> art_pt_set;
//    for(int i = 0; i < art_pt_vec.size(); ++i) {
//        art_pt_set.insert(art_pt_vec[i]);
//    }
//    
//    graph_hash_t::nodes_map_t::iterator it;
//    for(it = bcc.bcc_subgraph.nodes_map.begin();
//        it != bcc.bcc_subgraph.nodes_map.end();
//        ++it) {
//        node_id_t v = it->first;
//        if(art_pt_set.find(v) != art_pt_set.end()) {
//            // the node is an articulation point and it belongs to the bcc
//            // of interest, so it must be added to the art point of this bcc
//            // along with the sizes of the subgraphs it connects the bcc to
//            vector<int> subgraph_sz_vec;
//            vector<bool> visited_vec;
//            visited_vec.resize(size());
//            fill(visited_vec.begin(), visited_vec.end(), false);
//            visited_vec[v] = true;
//            vector<node_id_t> v_nbr_vec = get_nbrs(v);
//            for(int i = 0; i < v_nbr_vec.size(); ++i) {
//                node_id_t u = v_nbr_vec[i];
//                if(!bcc.bcc_subgraph.has_edge(v, u) && !visited_vec[u]) {
//                    //do a dfs from this guy to figure out the size
//                    //of the connected component connected to the bcc through v
//                    //and this edge
//                    int cnt = 0;
//                    stack<node_id_t> S;
//                    S.push(u);
//                    visited_vec[u] = true;
//                    while(!S.empty()) {
//                        node_id_t vv = S.top();
//                        S.pop();
//                        cnt++;
//                        vector<node_id_t> nbr_vec = get_nbrs(vv);
//                        for(int i = 0; i < nbr_vec.size(); ++i) {
//                            node_id_t nbr = nbr_vec[i];
//                            if(!visited_vec[nbr]) {
//                                visited_vec[nbr] = true;
//                                S.push(nbr);
//                            }
//                        }
//                    }
//                    subgraph_sz_vec.push_back(cnt);
//                }
//            }
//            bcc.art_pt_map.insert(make_pair(v, subgraph_sz_vec));
//        }
//    }
//}
//
//
//
//void bcc_scratch_t::compute_bc(tr1::unordered_map<node_id_t, double>& bc_map, int max_iter)
//{
//
//    //bc_map will have for each vertex in the muc it's betweenness centrality
//    //init bc of all nodes to zero
//    bcc_subgraph.i_fill_map<double>(bc_map, 0);
//    
//    //do BFS's from the nodes in the bcc
//    if(max_iter == -1) {
//        for(graph_hash_t::nodes_map_t::iterator
//                            it =  bcc_subgraph.nodes_map.begin();
//                            it != bcc_subgraph.nodes_map.end();
//                            ++it) {
//            node_id_t s = it->first;
//            bc_iter(s, bc_map);
//        }  
//    } else {
//        int cnt = 0;
//        for(graph_hash_t::nodes_map_t::iterator
//                            it =  bcc_subgraph.nodes_map.begin();
//                            it != bcc_subgraph.nodes_map.end();
//                            ++it) {
//            node_id_t s = it->first;
//            ++cnt;
//            bc_iter(s, bc_map);
//            if(cnt == max_iter)
//                break;
//        }  
//    }
//
//    if(max_iter == -1) {
//        for(art_pt_map_t::iterator
//            it =  art_pt_map.begin();
//            it != art_pt_map.end();
//            ++it) {
//            vector<int> size_vec;
//            size_vec = it->second;
//            if(size_vec.size() > 1) {
//                int sub = 0;
//                for(int i = 0; i < size_vec.size(); ++i) {
//                    sub += size_vec[i]*size_vec[i];
//                }
//                int VG_i = accumulate(it->second.begin(), it->second.end(), 0);
//                bc_map[it->first] = bc_map[it->first] + VG_i*VG_i - sub;
//            }
//        }
//        for (tr1::unordered_map<node_id_t, double>::iterator
//            it = bc_map.begin();
//                it != bc_map.end();
//                ++it) {
//            it->second = it->second / 2;
//        }    
//    }
//}
//
//
//void bcc_scratch_t::bc_iter(node_id_t s, tr1::unordered_map<node_id_t,double>& bc_map) {
//    tr1_map_t(vector<node_id_t>)&       P = bcc_subgraph.P;
//    tr1_map_t(int)&                     path_cnt_map = bcc_subgraph.path_cnt_vec;
//    tr1_map_t(int)&                     dist_map = bcc_subgraph.dist_vec;
//    tr1_map_t(double)&                  pair_dep_map = bcc_subgraph.pair_dep_vec;
//    tr1_map_t(double)&                  sigma_t_map = bcc_subgraph.sigma_t_map;
//    queue<node_id_t>&                   Q = bcc_subgraph.Q;
//    stack<node_id_t>                    S;
//
//    bcc_subgraph.init_maps();
//    path_cnt_map[s] = 1;
//    dist_map[s] = 0;
//
//    Q.push(s);
//    while (!Q.empty()) {
//        node_id_t v_i = Q.front();
//        Q.pop();
//        S.push(v_i);
//        for (int i = 0; i < bcc_subgraph.nodes_map[v_i].nbrs_vec.size(); ++i) {
//            node_id_t v_n = bcc_subgraph.nodes_map[v_i].nbrs_vec[i];
//            if (dist_map[v_n] < 0) {
//                Q.push(v_n);
//                dist_map[v_n] = dist_map[v_i] + 1;
//            }
//            if (dist_map[v_n] == dist_map[v_i] + 1) {
//                path_cnt_map[v_n] = path_cnt_map[v_n] + path_cnt_map[v_i];
//                P[v_n].push_back(v_i);
//            }
//        }
//    }
//    while (!S.empty()) {
//        node_id_t v_n = S.top();
//        S.pop();
//        if (art_pt_map.find(s) != art_pt_map.end()
//                && art_pt_map.find(v_n) != art_pt_map.end()
//                && s != v_n) {
//            int VG_s, VG_n;
//            VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
//            VG_n = accumulate(art_pt_map[v_n].begin(), art_pt_map[v_n].end(), 0);
//            int c_t = VG_s*VG_n;
//            sigma_t_map[v_n] = sigma_t_map[v_n] + c_t;
//            bc_map[v_n] = bc_map[v_n] + c_t;
//        }
//        for (int i = 0; i < P[v_n].size(); ++i) {
//            node_id_t v_p = P[v_n][i];
//            double sp_sn = ((double) path_cnt_map[v_p] / path_cnt_map[v_n]);
//            pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn * (1 + pair_dep_map[v_n]);
//            if (art_pt_map.find(s) != art_pt_map.end()) {
//                sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n] * sp_sn;
//                bc_map[v_p] = bc_map[v_p] + sigma_t_map[v_n] * sp_sn;
//            }
//        }
//        if (s != v_n) {
//            bc_map[v_n] = bc_map[v_n] + pair_dep_map[v_n];
//        }
//        if (art_pt_map.find(s) != art_pt_map.end()) {
//            int VG_s = accumulate(art_pt_map[s].begin(), art_pt_map[s].end(), 0);
//            bc_map[v_n] = bc_map[v_n] + pair_dep_map[v_n] * VG_s * 2;
//        }
//    }    
//}