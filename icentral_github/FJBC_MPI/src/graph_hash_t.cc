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

vector<node_id_t>& graph_hash_t::get_nbrs(node_id_t id)
{
    return nodes_map[id].nbrs_vec;
}


void graph_hash_t::init_maps()
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

template <typename T>
void graph_hash_t::i_fill_map(tr1_map_t(T)& map, T val)
{
    for(nodes_map_t::iterator it = nodes_map.begin();
            it != nodes_map.end();
            ++it) {
        node_id_t k = it->first;
        if(map.size() == 0)
            map.insert(make_pair(k, val));
        else
            map[k] = val;
    }
}

void graph_hash_t::insert_node(node_id_t id)
{
    if(nodes_map.find(id) == nodes_map.end()) {
        nodes_map.insert(make_pair(id, node_t()));
    }
}
void graph_hash_t::insert_edge(node_id_t src, node_id_t dst)
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
    
    insert_node(src);
    insert_node(dst);
    nodes_map[src].nbrs_vec.push_back(dst);
    nodes_map[dst].nbrs_vec.push_back(src);
}

void graph_hash_t::remove_edge(node_id_t src, node_id_t dst)
{
    pair<node_id_t, node_id_t> edge = make_pair(src, dst);
    edge_set.erase(edge);
    edge = make_pair(dst, src);
    edge_set.erase(edge);
    
    vector<node_id_t>& v1 = nodes_map[src].nbrs_vec;
    vector<node_id_t>& v2 = nodes_map[dst].nbrs_vec;
    v1.erase(std::remove(v1.begin(), v1.end(), dst), v1.end());
    v2.erase(std::remove(v2.begin(), v2.end(), src), v2.end());    
}

size_t graph_hash_t::size()
{
    return nodes_map.size();
}

void graph_hash_t::print_graph(bool with_nodes)
{
    printf("Edges: \n");
    for(set<edge_t >::iterator it = edge_set.begin(); it != edge_set.end(); it++) {
        printf("[%d] === [%d]\n", it->first, it->second);
    }
    if(with_nodes) {
        printf("Nodes with edge lists: \n");
        for(tr1::unordered_map<node_id_t, node_t>::iterator it = nodes_map.begin(); it != nodes_map.end(); it++) {
            printf("[%d]\t-- ", it->first);
            for(int i = 0; i < it->second.nbrs_vec.size(); i++) {
                printf("[%d] ", it->second.nbrs_vec[i]);
            }
            printf("\n");
        }
    }
}

bool graph_hash_t::has_edge(node_id_t src, node_id_t dst)
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
void graph_hash_t::conn_comp_sizes(vector<int>& out_vec)
{
    out_vec.clear();
    tr1::unordered_map<node_id_t, bool> visited_map;
    for(nodes_map_t::iterator it = nodes_map.begin(); it != nodes_map.end(); ++it) {
        visited_map.insert(make_pair(it->first, false));
    }
    for(tr1::unordered_map<node_id_t, bool>::iterator it = visited_map.begin(); it != visited_map.end(); ++it) {
        if(!it->second) {
            node_id_t s = it->first;
            it->second = true;
            out_vec.push_back(1);
            //do a BFS from s, count the number of vertices, and mark the visited
            queue<node_id_t> q;
            q.push(s);
            while(!q.empty()) {
                node_id_t node = q.front(); q.pop();
                for(int i = 0; i < nodes_map[node].nbrs_vec.size(); ++i) {
                    node_id_t nbr = nodes_map[node].nbrs_vec[i];
                    if(!visited_map[nbr]) {
                        visited_map[nbr] = true;
                        q.push(nbr);
                        out_vec.back() = out_vec.back()+1;
                    }
                }
            }
        }
    }
}

void graph_hash_t::find_conn_comp(vector<graph_hash_t>& out_vec)
{
    out_vec.clear();
    tr1_map_t(bool) visited_map;
    i_fill_map<bool>(visited_map, false);
    for(tr1::unordered_map<node_id_t, bool>::iterator it = visited_map.begin(); it != visited_map.end(); ++it) {
        if(!it->second) {
            graph_hash_t g;
            node_id_t s = it->first;
            it->second = true;
            g.insert_node(s);
            //do a BFS from s, count the number of vertices, and mark the visited
            queue<node_id_t> q;
            q.push(s);
            while(!q.empty()) {
                node_id_t node = q.front(); q.pop();
                for(int i = 0; i < nodes_map[node].nbrs_vec.size(); ++i) {
                    node_id_t nbr = nodes_map[node].nbrs_vec[i];
                    g.insert_edge(node, nbr);
                    if(!visited_map[nbr]) {
                        visited_map[nbr] = true;
                        q.push(nbr);
                    }
                }
            }
            out_vec.push_back(g);
        }
    }
}


void graph_hash_t::find_sssp(node_id_t s, tr1_map_t(int)& out_vec)
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


void graph_hash_t::fill_graph(graph_t& g)
{
    for(set<edge_t>::iterator it = g.edge_set.begin();
            it != g.edge_set.end();
            ++it) {
        insert_edge(it->first, it->second);
    }
    
}


void graph_hash_t::find_pruning_counts_exp(node_id_t src,
                                        node_id_t dst,
                                        int& d0_cnt,
                                        int& d1_cnt,
                                        int& d2_cnt)
{
    //IMP: the edge (src, dst) must not be in the graph
    //else doesn't make sense to count pruned BFS's
    bool insert = has_edge(src, dst);
    remove_edge(src, dst);
    tr1_map_t(int) d_src_vec, d_dst_vec;
    find_sssp(src, d_src_vec);
    find_sssp(dst, d_dst_vec);
    d0_cnt = 0;
    d1_cnt = 0;
    d2_cnt = 0;
    for(graph_hash_t::nodes_map_t::iterator
                        it =  nodes_map.begin();
                        it != nodes_map.end();
                        ++it) {
        node_id_t s = it->first;
        int diff = d_src_vec[s] - d_dst_vec[s];
        int abs_diff = abs(diff);
        switch(abs_diff) {
            case 0:  d0_cnt++; break;
            case 1:  d1_cnt++; break;
            default: d2_cnt++;
        }
    }
    if(insert)
        insert_edge(src, dst);
}