/* 
 * File:   graph_t.cc
 * Author: fuad
 * 
 * Created on Jun 23, 2014
 */
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include <stack>


#include "graph_t.h"
#include "mcb_find.h"
#include "bc.h"
#include "utility.h"

using namespace std;

vector<node_id_t>& subgraph_t::get_nbrs(node_id_t id)
{
    return nodes_vec[id];
}


void subgraph_t::init_maps()
{
    i_fill_map<vector<node_id_t> >(P, vector<node_id_t>());
    i_fill_map<int>(path_cnt_vec, 0);
    i_fill_map<int>(new_path_cnt_vec, 0);
    i_fill_map<int>(dist_vec, -1);
    i_fill_map<double>(pair_dep_vec, 0);
    i_fill_map<double>(new_pair_dep_vec, 0);
    i_fill_map<double>(sigma_t_map, 0);
    i_fill_map<double>(new_sigma_t_map, 0);
    i_fill_map<int>(path_cnt_inc_vec, 0);
    i_fill_map<bool>(visited_vec, false);
}


void subgraph_t::insert_edge(node_id_t src, node_id_t dst)
{
    edge_t e1, e2;
    e1.first   = src;
    e1.second  = dst;
    e2.first   = dst;
    e2.second  = src;
    if(edge_set.find(e1) != edge_set.end() || edge_set.find(e2) != edge_set.end())
        return;
    else
        edge_set.insert(e1);
    nodes_vec[src].push_back(dst);
    nodes_vec[dst].push_back(src);
}

void subgraph_t::remove_edge(node_id_t src, node_id_t dst)
{
    pair<node_id_t, node_id_t> edge = make_pair(src, dst);
    edge_set.erase(edge);
    edge = make_pair(dst, src);
    edge_set.erase(edge);
    
    vector<node_id_t>& v1 = nodes_vec[src];
    vector<node_id_t>& v2 = nodes_vec[dst];
    v1.erase(std::remove(v1.begin(), v1.end(), dst), v1.end());
    v2.erase(std::remove(v2.begin(), v2.end(), src), v2.end());    
}

size_t subgraph_t::size()
{
    return nodes_vec.size();
}

void subgraph_t::print_graph(bool with_nodes)
{
    printf("Edges: \n");
    for(set<edge_t >::iterator it = edge_set.begin(); it != edge_set.end(); it++) {
        printf("[%d] === [%d]\n", it->first, it->second);
    }
    if(with_nodes) {
        printf("Nodes with edge lists: \n");
        for(int i = 0; i < nodes_vec.size(); ++i) {
            printf("[%d]\t-- ", i);
            for(int j = 0; j < nodes_vec[i].size(); j++) {
                printf("[%d] ", nodes_vec[i][j]);
            }
            printf("\n");
        }
    }
}

bool subgraph_t::has_edge(node_id_t src, node_id_t dst)
{
    edge_t e1,e2;
    e1.first    = src;
    e1.second   = dst;
    e2.first    = dst;
    e2.second   = src;
    bool res = false;
    if(   edge_set.find(e1) != edge_set.end()
       || edge_set.find(e2) != edge_set.end()) {
        res = true;
    }
    return res;
}


//the vector @out_vec will have sizes of the connected
//components in the graph
void subgraph_t::conn_comp_sizes(vector<int>& out_vec)
{
    out_vec.clear();
    vector<bool> visited_vec;
    i_fill_map<bool>(visited_vec, false);
    for(int n_id = 0; n_id < visited_vec.size(); ++n_id) {
        if(!visited_vec[n_id]) {
            node_id_t s = n_id;
            visited_vec[n_id] = true;
            out_vec.push_back(1);
            //do a BFS from s, count the number of vertices, and mark the visited
            queue<node_id_t> q;
            q.push(s);
            while(!q.empty()) {
                node_id_t node = q.front(); q.pop();
                vector<node_id_t> nbr_vec = get_nbrs(node);
                for(int i = 0; i < nbr_vec.size(); ++i) {
                    node_id_t nbr = nbr_vec[i];
                    if(!visited_vec[nbr]) {
                        visited_vec[nbr] = true;
                        q.push(nbr);
                        out_vec.back() = out_vec.back()+1;
                    }
                }
            }
        }
    }
}


void subgraph_t::find_sssp(node_id_t s, vector<int>& out_vec)
{
    out_vec.clear();
    i_fill_map<int>(out_vec, -1);
    queue<node_id_t> Q;
    Q.push(s);
    out_vec[s] = 0;
    while(!Q.empty()) {
        node_id_t v = Q.front(); Q.pop();
        vector<node_id_t> nbr_vec = get_nbrs(v);
        for(int i = 0; i < nbr_vec.size(); ++i) {
            node_id_t nbr = nbr_vec[i];
            if(out_vec[nbr] == -1) {
                out_vec[nbr] = out_vec[v] + 1;
                Q.push(nbr);
            }
        }
    }
}


void subgraph_t::find_pruning_counts_exp(node_id_t src,
                                        node_id_t dst,
                                        int& d0_cnt,
                                        int& d1_cnt,
                                        int& d2_cnt)
{
    vector<int> d_src_vec, d_dst_vec;
    find_sssp(src, d_src_vec);
    find_sssp(dst, d_dst_vec);
    d0_cnt = 0;
    d1_cnt = 0;
    d2_cnt = 0;
    for(int i = 0; i < nodes_vec.size(); ++i) {
        node_id_t s = i;
        int diff = d_src_vec[s] - d_dst_vec[s];
        int abs_diff = abs(diff);
        switch(abs_diff) {
            case 0:  d0_cnt++; break;
            case 1:  d1_cnt++; break;
            default: d2_cnt++;
        }
    }
}


void subgraph_t::fill_graph(graph_hash_t& g)
{       
    nodes_vec.clear();
    edge_set.clear();
    
    inout_label_map.clear(); // maps internal labels to original
    outin_label_map.clear(); // original labels to internal    

    nodes_vec.resize(g.size());
    inout_label_map.resize(g.size());
    int cnt = 0;
    for(graph_hash_t::nodes_map_t::iterator it = g.nodes_map.begin();
            it != g.nodes_map.end();
            ++it) {
        inout_label_map[cnt] = it->first;
        outin_label_map.insert(make_pair(it->first, cnt));
        cnt++;
    }
    for(set<edge_t>::iterator it = g.edge_set.begin();
            it != g.edge_set.end();
            ++it) {
        node_id_t src, dst;
        src = outin_label_map[it->first];
        dst = outin_label_map[it->second];
        insert_edge(src, dst);
    }
}

void subgraph_t::fill_graph(graph_t& g)
{       
    nodes_vec.clear();
    edge_set.clear();
    
    inout_label_map.clear(); // maps internal labels to original
    outin_label_map.clear(); // original labels to internal    

    nodes_vec.resize(g.size());
    for(node_id_t i = 0; i < g.size(); ++i) {
        nodes_vec[i] = g.nodes_vec[i].nbrs_vec;
    } 
    edge_set = g.edge_set;
    
    inout_label_map.resize(size());
    for(node_id_t i = 0; i < size(); ++i) {
        outin_label_map.insert(make_pair(i, i));
    }
}
