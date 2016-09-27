#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <stack>
#include <queue>

#include "bc.h"
#include "utility.h"

using namespace std;

vector<double> naive_betweenness_centrality(graph_t& graph_in)
{
    vector<double> bc_vec;
    bc_vec.resize(graph_in.size());
    
    /*1. find all-pair shortest paths and store in a n^2 array*/
    /*2. find all-pair shortest paths count and store in a x^2 array*/
    /*3. do the three nested loops to find BC of each node*/
    
    vector<vector<int> > dst_2vec;      //all-pair distance array
    vector<vector<int> > cnt_2vec;      //all-pair shortest path count array
    dst_2vec.resize(graph_in.size());
    cnt_2vec.resize(graph_in.size());
    
    vector<int> tmp_vec;
    tmp_vec.resize(graph_in.size());
    
 
    /*Floyd-Warshall to compute all-pair shortest paths*/
    /*init distance with (infinity)*/
    fill(tmp_vec.begin(), tmp_vec.end(), INF);
    fill(dst_2vec.begin(), dst_2vec.end(), tmp_vec);
    /*update distance for immediate neighbors*/
    for(node_id_t i = 0; i < graph_in.size(); i++) {
        dst_2vec[i][i] = 0;
        for(int j = 0; j < graph_in.get_nbrs(i).size(); j++) {
            node_id_t nbr_id = graph_in.get_nbrs(i)[j];
            dst_2vec[i][nbr_id] = 1;
        }
    }
    for(int k = 0; k < graph_in.size(); ++k) {
        for(int i = 0; i < dst_2vec.size(); ++i) {
            for(int j = 0; j < dst_2vec.size(); ++j) {
                if(dst_2vec[i][j] > dst_2vec[i][k] + dst_2vec[k][j]) {
                    dst_2vec[i][j] = dst_2vec[i][k] + dst_2vec[k][j];
                }
            }
        }
    }
    
    /*Floyd-Warshall like algorithm to compute all-pair shortest paths counts*/
    /*init counts with 0*/
    fill(tmp_vec.begin(), tmp_vec.end(), 0);
    fill(cnt_2vec.begin(), cnt_2vec.end(), tmp_vec);
    /*update counts for immediate neighbors*/
    for(node_id_t i = 0; i < graph_in.size(); i++) {
        cnt_2vec[i][i] = 0;
        for(int j = 0; j < graph_in.get_nbrs(i).size(); j++) {
            node_id_t nbr_id = graph_in.get_nbrs(i)[j];
            cnt_2vec[i][nbr_id] = 1;
        }
    }
    for(int k = 0; k < graph_in.size(); ++k) {
        for(int i = 0; i < cnt_2vec.size(); ++i) {
            for(int j = 0; j < cnt_2vec.size(); ++j) {
                if(dst_2vec[i][j] == dst_2vec[i][k] + dst_2vec[k][j]) {
                    cnt_2vec[i][j] += cnt_2vec[i][k] * cnt_2vec[k][j];
                }
            }
        }
    }
    /*betweenness centrality computation*/
    for(node_id_t s = 0; s < graph_in.size(); ++s) {
        for(node_id_t t = 0; t < graph_in.size(); ++t) {
            for(node_id_t v = 0; v < bc_vec.size(); ++v) {
                if(s != v && v != t && s != t) {
                    if(dst_2vec[s][t] == dst_2vec[s][v] + dst_2vec[v][t]) {
                        bc_vec[v] += (cnt_2vec[s][v]*cnt_2vec[v][t]/(double)cnt_2vec[s][t]);
                    }
                }
            }
        }
    }
    
    for(node_id_t i = 0; i < bc_vec.size(); ++i) {
        bc_vec[i] = bc_vec[i]/2;
    }     
    return bc_vec;
}

vector<double> brandes_betweenness_centrality(graph_t& graph_in)
{
    vector<double> bc_vec;
    bc_vec.resize(graph_in.size());
    fill(bc_vec.begin(), bc_vec.end(), 0.0);
    
    for(node_id_t s = 0; s < graph_in.size(); ++s) {
        stack<node_id_t>                  S;
        vector<vector<node_id_t> >        P;
        vector<int>                     path_cnt_vec;
        vector<int>                     dist_vec;
        queue<node_id_t>                  Q;
        
        P.resize(graph_in.size());
        path_cnt_vec.resize(graph_in.size());
        fill(path_cnt_vec.begin(), path_cnt_vec.end(), 0);
        path_cnt_vec[s] = 1;
        dist_vec.resize(graph_in.size());
        fill(dist_vec.begin(), dist_vec.end(), -1);
        dist_vec[s] = 0;
        
        Q.push(s);
        
        while(!Q.empty()) {
            node_id_t v = Q.front(); Q.pop();
            S.push(v);
            vector<node_id_t> nbr_vec = graph_in.get_nbrs(v);
            for(int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
                node_id_t w = nbr_vec[nbr_idx];
                if(dist_vec[w] < 0) {
                    Q.push(w);
                    dist_vec[w] = dist_vec[v] + 1;
                }
                if(dist_vec[w] == dist_vec[v] + 1) {
                    path_cnt_vec[w] += path_cnt_vec[v];
                    P[w].push_back(v);
                }
            }
        }
        vector<double> pair_dep_vec;
        pair_dep_vec.resize(graph_in.size());
        fill(pair_dep_vec.begin(), pair_dep_vec.end(), 0.0);
        while(!S.empty()) {
            node_id_t w = S.top(); S.pop();
            for(int idx = 0; idx < P[w].size(); ++idx) {
                node_id_t v = P[w][idx];
                pair_dep_vec[v] += (((double)path_cnt_vec[v]/path_cnt_vec[w]) * (1+pair_dep_vec[w]));
            }
            if (w != s) {
                bc_vec[w] += pair_dep_vec[w];
            }
        }
    }
    for(node_id_t i = 0; i < bc_vec.size(); ++i) {
        bc_vec[i] = bc_vec[i]/2;
    }
    return bc_vec;
}

vector<double> brandes_bc(graph_t& graph_in, int max_iter)
{
    vector<double> bc_vec;
    bc_vec.resize(graph_in.size());
    fill(bc_vec.begin(), bc_vec.end(), 0.0);
    if(max_iter == -1) {
        for(node_id_t s = 0; s < graph_in.size(); ++s) {
            brandes_iter(graph_in, s, bc_vec);
        }
        for(node_id_t i = 0; i < bc_vec.size(); ++i) {
            bc_vec[i] = bc_vec[i]/2;
        }
    } else {
        int num_iter = min(max_iter, (int)graph_in.size());
        for(node_id_t s = 0; s < num_iter; ++s) {
            brandes_iter(graph_in, s, bc_vec);
        }
    }

    return bc_vec;
}

//@s_out_vec must have the same size as the graph
void brandes_iter(graph_t& graph_in,
        node_id_t s,
        vector<double>& bc_vec)
{
    stack<node_id_t> S;
    vector<vector<node_id_t> > P;
    vector<int> path_cnt_vec;
    vector<int> dist_vec;
    queue<node_id_t> Q;
    vector<double> pair_dep_vec;

    P.resize(graph_in.size());
    path_cnt_vec.resize(graph_in.size());
    fill(path_cnt_vec.begin(), path_cnt_vec.end(), 0);
    path_cnt_vec[s] = 1;
    dist_vec.resize(graph_in.size());
    fill(dist_vec.begin(), dist_vec.end(), -1);
    dist_vec[s] = 0;
    pair_dep_vec.resize(graph_in.size());
    fill(pair_dep_vec.begin(), pair_dep_vec.end(), 0.0);

    Q.push(s);

    while (!Q.empty()) {
        node_id_t v = Q.front();
        Q.pop();
        S.push(v);
        vector<node_id_t> nbr_vec = graph_in.get_nbrs(v);
        for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
            node_id_t w = nbr_vec[nbr_idx];
            if (dist_vec[w] < 0) {
                Q.push(w);
                dist_vec[w] = dist_vec[v] + 1;
            }
            if (dist_vec[w] == dist_vec[v] + 1) {
                path_cnt_vec[w] += path_cnt_vec[v];
                P[w].push_back(v);
            }
        }
    }

    while (!S.empty()) {
        node_id_t w = S.top();
        S.pop();
        for (int idx = 0; idx < P[w].size(); ++idx) {
            node_id_t v = P[w][idx];
            pair_dep_vec[v] += (((double) path_cnt_vec[v] / path_cnt_vec[w]) * (1 + pair_dep_vec[w]));
        }
        if (w != s) {
            bc_vec[w] += pair_dep_vec[w];
        }
    }
}


tr1_map_t(double) brandes_bc_hash_out(graph_t& graph_in, int max_iter)
{
    tr1_map_t(double) bc_vec;
    for(int i = 0; i < graph_in.nodes_vec.size(); ++i) {
        bc_vec.insert(make_pair(i, 0.0));
    }

    if(max_iter == -1) {
        for(node_id_t s = 0; s < graph_in.size(); ++s) {
            brandes_iter_hash_out(graph_in, s, bc_vec);
        }
        for(node_id_t i = 0; i < bc_vec.size(); ++i) {
            bc_vec[i] = bc_vec[i]/2;
        }
    } else {
        int num_iter = min(max_iter, (int)graph_in.size());
        for(node_id_t s = 0; s < num_iter; ++s) {
            brandes_iter_hash_out(graph_in, s, bc_vec);
        }
    }

    return bc_vec;
}

void brandes_iter_hash_out(graph_t& graph_in,
        node_id_t s,
        tr1_map_t(double)& bc_vec)
{
    stack<node_id_t> S;
    vector<vector<node_id_t> > P;
    vector<int> path_cnt_vec;
    vector<int> dist_vec;
    queue<node_id_t> Q;
    vector<double> pair_dep_vec;

    P.resize(graph_in.size());
    path_cnt_vec.resize(graph_in.size());
    fill(path_cnt_vec.begin(), path_cnt_vec.end(), 0);
    path_cnt_vec[s] = 1;
    dist_vec.resize(graph_in.size());
    fill(dist_vec.begin(), dist_vec.end(), -1);
    dist_vec[s] = 0;
    pair_dep_vec.resize(graph_in.size());
    fill(pair_dep_vec.begin(), pair_dep_vec.end(), 0.0);

    Q.push(s);

    while (!Q.empty()) {
        node_id_t v = Q.front();
        Q.pop();
        S.push(v);
        vector<node_id_t> nbr_vec = graph_in.get_nbrs(v);
        for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
            node_id_t w = nbr_vec[nbr_idx];
            if (dist_vec[w] < 0) {
                Q.push(w);
                dist_vec[w] = dist_vec[v] + 1;
            }
            if (dist_vec[w] == dist_vec[v] + 1) {
                path_cnt_vec[w] += path_cnt_vec[v];
                P[w].push_back(v);
            }
        }
    }

    while (!S.empty()) {
        node_id_t w = S.top();
        S.pop();
        for (int idx = 0; idx < P[w].size(); ++idx) {
            node_id_t v = P[w][idx];
            pair_dep_vec[v] += (((double) path_cnt_vec[v] / path_cnt_vec[w]) * (1 + pair_dep_vec[w]));
        }
        if (w != s) {
            bc_vec[w] += pair_dep_vec[w];
        }
    }
}

void incremental_brandes(graph_t& graph_in,
        node_id_t src,
        node_id_t dst,
        vector<double>& bc_vec)
{
    vector<int> d_src_vec, d_dst_vec;
    d_src_vec.resize(graph_in.size());
    d_dst_vec.resize(graph_in.size());
    fill(d_src_vec.begin(), d_src_vec.end(), -1);
    fill(d_dst_vec.begin(), d_dst_vec.end(), -1);
    graph_in.find_sssp(src, d_src_vec);
    graph_in.find_sssp(dst, d_dst_vec);
    
    //fill(bc_vec.begin(), bc_vec.end(), 0);
    
    for(node_id_t s = 0; s < graph_in.size(); ++s) {
        if(d_src_vec[s] != d_dst_vec[s]) { 
            brandes_delta_iter(graph_in,
                    s, 
                    src, 
                    dst, 
                    d_src_vec[s], 
                    d_dst_vec[s], 
                    bc_vec);
        }
    }
}

void incremental_brandes_experimental(graph_t& graph_in,
        node_id_t src,
        node_id_t dst,
        vector<double>& bc_vec,
        vector<double>& time_vec,
        vector<int>&    cnt_vec)
{
    vector<int> d_src_vec, d_dst_vec;
    d_src_vec.resize(graph_in.size());
    d_dst_vec.resize(graph_in.size());
    fill(d_src_vec.begin(), d_src_vec.end(), -1);
    fill(d_dst_vec.begin(), d_dst_vec.end(), -1);
    graph_in.find_sssp(src, d_src_vec);
    graph_in.find_sssp(dst, d_dst_vec);
    
    //fill(bc_vec.begin(), bc_vec.end(), 0);
    
    int    cnt_arr[] = {0,0,0};
    double tot_arr[] = {0,0,0};
    for(node_id_t s = 0; s < graph_in.size(); ++s) {
        timer tm;
        tm.start();
        if(d_src_vec[s] != d_dst_vec[s]) { 
            brandes_delta_iter(graph_in,
                    s, 
                    src, 
                    dst, 
                    d_src_vec[s], 
                    d_dst_vec[s], 
                    bc_vec);
        }
        tm.stop();
        int diff = d_src_vec[s] - d_dst_vec[s];
        int abs_diff = abs(diff);
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
    }
    cnt_vec.insert(cnt_vec.begin(), cnt_arr, cnt_arr+3);
    time_vec.insert(time_vec.begin(), tot_arr, tot_arr+3);
}

//the two cases, when d(s, dst) - d(s, src) is 1 and > 2 are handled here
void brandes_delta_iter(graph_t& graph_in,
        node_id_t s,
        node_id_t src,
        node_id_t dst,
        int d_src,
        int d_dst,
        vector<double>& bc_vec)
{
    //make sure that @src is the closer to source node
    if(d_src > d_dst) {
        node_id_t       tmp_v;
        int             tmp_i;
        tmp_v = src; src = dst; dst = tmp_v;
        tmp_i = d_src; d_src = d_dst; d_dst = tmp_i;
        
    }
    
    //case 1: d_dst - d_src = 1 (the sort of easier case)
    if(d_dst - d_src == 1) {//if(d_dst - d_src == 1) {
        stack<node_id_t> S;
        vector<vector<node_id_t> > P;
        vector<int> old_path_cnt_vec;
        vector<int> new_path_cnt_vec; //NEW
        vector<int> dist_vec;
        queue<node_id_t> Q;
        vector<double> old_pair_dep_vec;
        vector<double> new_pair_dep_vec;//NEW

        P.resize(graph_in.size());
        new_path_cnt_vec.resize(graph_in.size());
        old_path_cnt_vec.resize(graph_in.size());
        dist_vec.resize(graph_in.size());
        fill(old_path_cnt_vec.begin(), old_path_cnt_vec.end(), 0);
        fill(new_path_cnt_vec.begin(), new_path_cnt_vec.end(), 0);
        fill(dist_vec.begin(), dist_vec.end(), -1);
        old_path_cnt_vec[s] = 1;
        new_path_cnt_vec[s] = 1;
        dist_vec[s] = 0;
        old_pair_dep_vec.resize(graph_in.size());
        new_pair_dep_vec.resize(graph_in.size());
        fill(old_pair_dep_vec.begin(), old_pair_dep_vec.end(), 0.0);
        fill(new_pair_dep_vec.begin(), new_pair_dep_vec.end(), 0.0);

        Q.push(s);

        while (!Q.empty()) {
            node_id_t v = Q.front();
            Q.pop();
            S.push(v);
            if(v == dst) {//NEW
                new_path_cnt_vec[v] += old_path_cnt_vec[src];
            }
            vector<node_id_t> nbr_vec = graph_in.get_nbrs(v);
            for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
                node_id_t w = nbr_vec[nbr_idx];
                if (dist_vec[w] < 0) {
                    Q.push(w);
                    dist_vec[w] = dist_vec[v] + 1;
                }
                if (dist_vec[w] == dist_vec[v] + 1) {
                    old_path_cnt_vec[w] += old_path_cnt_vec[v];
                    new_path_cnt_vec[w] += new_path_cnt_vec[v];//NEW
                    P[w].push_back(v);
                }
            }
        }

        while (!S.empty()) {
            node_id_t w = S.top();
            S.pop();
            for (int idx = 0; idx < P[w].size(); ++idx) {
                node_id_t v = P[w][idx];
                old_pair_dep_vec[v] += (((double) old_path_cnt_vec[v] / old_path_cnt_vec[w]) * (1 + old_pair_dep_vec[w]));
                new_pair_dep_vec[v] += (((double) new_path_cnt_vec[v] / new_path_cnt_vec[w]) * (1 + new_pair_dep_vec[w]));
            }
            //IMP: this is the only change that happens to P, @src should be added as parent for dst
            if(w == dst) {
                node_id_t v = src;
                new_pair_dep_vec[v] += (((double) new_path_cnt_vec[v] / new_path_cnt_vec[w]) * (1 + new_pair_dep_vec[w]));
            }
            if (w != s) {
                bc_vec[w] += (new_pair_dep_vec[w] - old_pair_dep_vec[w])/2.0;
                //bc_vec[w] += (0 - old_pair_dep_vec[w])/2.0;
            }
        }
    } else {
        vector<vector<node_id_t> > P;
        vector<int> path_cnt_vec;
        vector<int> path_cnt_inc_vec;
        vector<int> dist_vec;
        vector<double> pair_dep_vec;
        vector<bool>& visited_vec = graph_in.tmp_bool_vec;
        queue<node_id_t> Q;
        vector<node_id_t> S;
        
        P.resize(graph_in.size());
        dist_vec.resize(graph_in.size());
        path_cnt_vec.resize(graph_in.size());
        pair_dep_vec.resize(graph_in.size());
        path_cnt_inc_vec.resize(graph_in.size());
        
        fill(path_cnt_vec.begin(), path_cnt_vec.end(), 0);
        path_cnt_vec[s] = 1;
        fill(dist_vec.begin(), dist_vec.end(), -1);
        dist_vec[s] = 0;
        fill(pair_dep_vec.begin(), pair_dep_vec.end(), 0.0);
        fill(path_cnt_inc_vec.begin(), path_cnt_inc_vec.end(), 0);
        
        Q.push(s);
        while (!Q.empty()) {
            node_id_t v = Q.front();
            Q.pop();
            //S.push(v);
            S.push_back(v);
            vector<node_id_t> nbr_vec = graph_in.get_nbrs(v);
            for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
                node_id_t w = nbr_vec[nbr_idx];
                if (dist_vec[w] < 0) {
                    Q.push(w);
                    dist_vec[w] = dist_vec[v] + 1;
                }
                if (dist_vec[w] == dist_vec[v] + 1) {
                    path_cnt_vec[w] += path_cnt_vec[v];
                    P[w].push_back(v);
                }
            }
        }

        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t w = S[i];
            for (int idx = 0; idx < P[w].size(); ++idx) {
                node_id_t v = P[w][idx];
                pair_dep_vec[v] += (((double) path_cnt_vec[v] / path_cnt_vec[w]) * (1 + pair_dep_vec[w]));
            }
            if (w != s) {
                bc_vec[w] -= pair_dep_vec[w]/2.0;
            }
        }

        //vector<int> old_dist_vec = dist_vec;    
        Q.push(dst);
        dist_vec[dst] = dist_vec[src] + 1;
        P[dst].clear();
        P[dst].push_back(src);
        path_cnt_inc_vec[dst] = path_cnt_vec[src];
        path_cnt_vec[dst] = path_cnt_inc_vec[dst];
        //tr1::unordered_set<node_id_t> visited_set;
        
        //visited_vec.resize(graph_in.size());
        fill(visited_vec.begin(), visited_vec.end(), false);
        visited_vec[dst] = true;
        //visited_set.insert(dst);
        
        while (!Q.empty()) {
            node_id_t v = Q.front();
            Q.pop();
            //S.push(v);
            vector<node_id_t> nbr_vec = graph_in.get_nbrs(v);
            for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
                node_id_t w = nbr_vec[nbr_idx];
                if (dist_vec[w] > (dist_vec[v] + 1)) {
                    dist_vec[w] = dist_vec[v] + 1;
                    P[w].clear();
                    P[w].push_back(v);
                    path_cnt_vec[w] = 0;
                    path_cnt_inc_vec[w] = path_cnt_inc_vec[v];
                    path_cnt_vec[w] += path_cnt_inc_vec[w];
//                    if(!visited_vec[w]) {
//                        visited_vec[w] = true;
//                        Q.push(w);
//                    }
                } else if (dist_vec[w] == (dist_vec[v] + 1)) {
                    path_cnt_inc_vec[w] += path_cnt_inc_vec[v];
                    path_cnt_vec[w] += path_cnt_inc_vec[v];
                    //if(old_dist_vec[w] == old_dist_vec[v] || v == dst) {
                    if(find(P[w].begin(), P[w].end(), v) == P[w].end()) {
                        P[w].push_back(v);
                    } 
//                    if(!visited_vec[w]) {
//                        visited_vec[w] = true;
//                        Q.push(w);
//                    }
                }
                if(!visited_vec[w]) {
                    visited_vec[w] = true;
                    Q.push(w);
                }
            }
        }
        
        
        //IMP::THIS CAN BE MADE MUCH BETTER!
        //HEAP FOR EXAMPLE
        //EVEN THE SWAPPING CAN BE DONE MORE EFFICIENTLY
        for(int i = 1; i < S.size(); ++i) {
            if(dist_vec[S[i-1]] > dist_vec[S[i]]) {
                int j = i;
                while(dist_vec[S[j-1]] > dist_vec[S[j]]) {
                    node_id_t tmp = S[j-1];
                    S[j-1] = S[j];
                    S[j] = tmp;
                    --j;
                }
            }
        }
        
        fill(pair_dep_vec.begin(), pair_dep_vec.end(), 0.0);
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t w = S[i];
            for (int idx = 0; idx < P[w].size(); ++idx) {
                node_id_t v = P[w][idx];
                pair_dep_vec[v] += (((double) path_cnt_vec[v] / path_cnt_vec[w]) * (1 + pair_dep_vec[w]));
            }
            if (w != s) {
                bc_vec[w] += pair_dep_vec[w]/2.0;
            }
        }
    }
}

tr1_map_t(double) brandes_bc_hash(graph_hash_t& graph_in, int max_iter)
{
    tr1_map_t(double) bc_vec;
    graph_in.i_fill_map<double>(bc_vec, 0);
    //means do all iterations
    if(max_iter == -1) {
        for(graph_hash_t::nodes_map_t::iterator it = graph_in.nodes_map.begin();
                it != graph_in.nodes_map.end(); 
                ++it) {
            node_id_t s = it->first;
            brandes_iter(graph_in, s, bc_vec);
        }
    } else {
        int cnt = 0;
        for(graph_hash_t::nodes_map_t::iterator it = graph_in.nodes_map.begin();
                it != graph_in.nodes_map.end(); 
                ++it) {
            node_id_t s = it->first;
            ++cnt;
            brandes_iter(graph_in, s, bc_vec);
            if(cnt == max_iter)
                break;
        }
    }
    
    if(max_iter == -1) {
        for(tr1_map_t(double)::iterator it = bc_vec.begin();
                it != bc_vec.end();
                ++it) {
            it->second = it->second/2;
        }
    }
    return bc_vec;    
}

void brandes_iter(graph_hash_t& graph_in,
        node_id_t s,
        tr1_map_t(double)& bc_vec)
{
    tr1_map_t(vector<node_id_t> )&      P = graph_in.P;
    tr1_map_t(int)&                     path_cnt_vec = graph_in.path_cnt_vec;
    tr1_map_t(int)&                     dist_vec = graph_in.dist_vec;
    tr1_map_t(double)&                  pair_dep_vec = graph_in.pair_dep_vec;
    vector<node_id_t>&                  S = graph_in.S;
    queue<node_id_t>&                   Q = graph_in.Q;
    
    graph_in.init_maps();
    
    S.clear();
    
    path_cnt_vec[s] = 1;
    dist_vec[s] = 0;

    Q.push(s);

    while (!Q.empty()) {
        node_id_t v = Q.front();
        Q.pop();
        S.push_back(v);
        vector<node_id_t> nbr_vec = graph_in.get_nbrs(v);
        for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
            node_id_t w = nbr_vec[nbr_idx];
            if (dist_vec[w] < 0) {
                Q.push(w);
                dist_vec[w] = dist_vec[v] + 1;
            }
            if (dist_vec[w] == dist_vec[v] + 1) {
                path_cnt_vec[w] += path_cnt_vec[v];
                P[w].push_back(v);
            }
        }
    }

    for(int i = S.size()-1; i >= 0; --i) {
        node_id_t w = S[i];
        for(int idx = 0; idx < P[w].size(); ++idx) {
            node_id_t v = P[w][idx];
            pair_dep_vec[v] += (((double) path_cnt_vec[v] / path_cnt_vec[w]) * (1 + pair_dep_vec[w]));
        }
        if (w != s) {
            bc_vec[w] += pair_dep_vec[w];
        }
    }    
}