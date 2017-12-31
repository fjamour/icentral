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


#include "graph_t.h"

/*
 * This function will return a vector of the articulation points,
 * if an articulation point appears in more than two biconnected
 * components it will appear more than once (number of bcc's it appears in - 1)
 */
void graph_t::find_art_points(vector<node_id_t>& out_vec)
{
    node_id_t u;
    vector<int> color_vec;
    vector<int> low_vec;
    vector<int> d_vec;
    vector<node_id_t> pred_vec;
    int time;
    
    u = 0;
    color_vec.resize(size());
    low_vec.resize(size());
    d_vec.resize(size());
    pred_vec.resize(size());
    fill(color_vec.begin(), color_vec.end(), 0);
    fill(low_vec.begin(), low_vec.end(), -1);
    fill(d_vec.begin(), d_vec.end(), -1);
    fill(pred_vec.begin(), pred_vec.end(), -1);
    
    i_art_point_dfs_visitor(u, color_vec, low_vec, d_vec, pred_vec, time, out_vec);
}

void graph_t::i_art_point_dfs_visitor(node_id_t u,
                          vector<int>& color_vec,
                          vector<int>& low_vec,
                          vector<int>& d_vec,
                          vector<node_id_t>& pred_vec,
                          int& time,
                          vector<node_id_t>& art_point_vec)
{
    color_vec[u] = 1;           //1 means grey
    d_vec[u] = ++time;
    low_vec[u] = d_vec[u];
    vector<node_id_t> nbr_vec = get_nbrs(u);
    int tree_edge_cnt = 0;
    for(int i = 0; i < nbr_vec.size(); ++i) {
        node_id_t v = nbr_vec[i];
        if(color_vec[v] == 0) { // (u, v) is a tree edge
            pred_vec[v] = u;
            tree_edge_cnt++;
            i_art_point_dfs_visitor(v, color_vec, low_vec, d_vec, pred_vec, time, art_point_vec);
            low_vec[u] = min(low_vec[u], low_vec[v]);
            if(pred_vec[u] == -1) {     //this is the root of the tree
                //if v is u's second child then u is art vertex
                if(tree_edge_cnt > 1) {
                    art_point_vec.push_back(u);
                }
            } else if(low_vec[v] >= d_vec[u]) {
                art_point_vec.push_back(u);
            }
        } else if(v != pred_vec[u]) {   // (u, v) is a back edge
            low_vec[u] = min(low_vec[u], d_vec[v]);
        }
    }
}

void graph_t::find_bicon_comp(vector<graph_hash_t>& out_vec)
{
    node_id_t u;
    vector<int> color_vec;
    vector<int> low_vec;
    vector<int> d_vec;
    vector<node_id_t> pred_vec;
    stack<edge_t> edge_stack;
    int time;
    
    u = 0;
    color_vec.resize(size());
    low_vec.resize(size());
    d_vec.resize(size());
    pred_vec.resize(size());
    fill(color_vec.begin(), color_vec.end(), 0);
    fill(low_vec.begin(), low_vec.end(), -1);
    fill(d_vec.begin(), d_vec.end(), -1);
    fill(pred_vec.begin(), pred_vec.end(), -1);
    
    i_bcc_dfs_visitor(u,
                    color_vec,
                    low_vec,
                    d_vec,
                    pred_vec,
                    edge_stack,
                    time,
                    out_vec);   
}

void graph_t::i_bcc_dfs_visitor(node_id_t u,
                                    vector<int>& color_vec,
                                    vector<int>& low_vec,
                                    vector<int>& d_vec,
                                    vector<node_id_t>& pred_vec,
                                    stack<edge_t>& edge_stack,
                                    int& time,
                                    vector<graph_hash_t>& bcc_vec)
{
    color_vec[u] = 1;           //1 means grey
    d_vec[u] = ++time;
    low_vec[u] = d_vec[u];
    vector<node_id_t> nbr_vec = get_nbrs(u);
    int tree_edge_cnt = 0;
    for(int i = 0; i < nbr_vec.size(); ++i) {
        node_id_t v = nbr_vec[i];
        if(color_vec[v] == 0) { // (u, v) is a tree edge
            edge_stack.push(make_pair(u, v));
            pred_vec[v] = u;
            tree_edge_cnt++;
            i_bcc_dfs_visitor(v,
                    color_vec,
                    low_vec,
                    d_vec,
                    pred_vec,
                    edge_stack,
                    time,
                    bcc_vec);
            low_vec[u] = min(low_vec[u], low_vec[v]);
            if(low_vec[v] >= d_vec[u]) {
                // u is an articulation point, ready to output a bcc
                edge_t last_edge = make_pair(u, v);
                edge_t e;
                graph_hash_t bcc;
                do {
                    e = edge_stack.top();
                    edge_stack.pop();
                    bcc.insert_edge(e.first, e.second);
                } while(e != last_edge);
                bcc_vec.push_back(bcc);
            }
        } else if(v != pred_vec[u] && d_vec[v] < d_vec[u]) {   // (u, v) is a back edge
            edge_stack.push(make_pair(u, v));
            low_vec[u] = min(low_vec[u], d_vec[v]);
        }
    }
}


// the edge (@src, @dst) must be in the graph
// the biconnected component subgraph that contains the passed edge
// will be returned in @bcc
// IMP assumption:
// the input graph is connected, and the edge (src, dst) is inserted,
// so the edge (src, dst) is part of a cycle, so both ends belong to the
// same bcc
//IMP: (src,dst) must exist!
void graph_t::find_edge_bcc_subgraph(graph_hash_t& bcc, node_id_t src, node_id_t dst)
{
    node_id_t u;
    vector<int> color_vec;
    vector<int> low_vec;
    vector<int> d_vec;
    vector<node_id_t> pred_vec;
    stack<edge_t> edge_stack;
    int time;
    
    color_vec.resize(size());
    low_vec.resize(size());
    d_vec.resize(size());
    pred_vec.resize(size());
    fill(color_vec.begin(), color_vec.end(), 0);
    fill(low_vec.begin(), low_vec.end(), -1);
    fill(d_vec.begin(), d_vec.end(), -1);
    fill(pred_vec.begin(), pred_vec.end(), -1);
    
    //(src, dst) needs to be the first edge to recurse from
    //IMP: explain this clearly
    //Jun 16, 2014.. I have no idea why!
    for(int i = 0; i < nodes_vec[src].nbrs_vec.size(); ++i) {
        if(nodes_vec[src].nbrs_vec[i] == dst) {
            nodes_vec[src].nbrs_vec[i] = nodes_vec[src].nbrs_vec[0];
            nodes_vec[src].nbrs_vec[0] = dst;
        }
    }
    
    i_edge_bcc_dfs_visitor(src,
                    color_vec,
                    low_vec,
                    d_vec,
                    pred_vec,
                    edge_stack,
                    time);
      
    while(!edge_stack.empty()) {
        edge_t e = edge_stack.top();
        edge_stack.pop();
        bcc.insert_edge(e.first, e.second);
    }
}

void graph_t::i_edge_bcc_dfs_visitor(node_id_t u,
                                    vector<int>& color_vec,
                                    vector<int>& low_vec,
                                    vector<int>& d_vec,
                                    vector<node_id_t>& pred_vec,
                                    stack<edge_t>& edge_stack,
                                    int& time)
{
    color_vec[u] = 1;           //1 means grey
    d_vec[u] = ++time;
    low_vec[u] = d_vec[u];
    vector<node_id_t> nbr_vec = get_nbrs(u);
    int tree_edge_cnt = 0;
    for(int i = 0; i < nbr_vec.size(); ++i) {
        node_id_t v = nbr_vec[i];
        if(color_vec[v] == 0) { // (u, v) is a tree edge
            edge_stack.push(make_pair(u, v));
            pred_vec[v] = u;
            tree_edge_cnt++;
            /// f() starts
            i_edge_bcc_dfs_visitor(v,
                    color_vec,
                    low_vec,
                    d_vec,
                    pred_vec,
                    edge_stack,
                    time);
            /// g() starts
            low_vec[u] = min(low_vec[u], low_vec[v]);
            bool flush_edges = (low_vec[v] >= d_vec[u] && pred_vec[u] != -1);
            flush_edges |= (pred_vec[u] == -1 && tree_edge_cnt >= 2);
            if(flush_edges) {
                // u is an articulation point, ready to output a bcc
                edge_t last_edge = make_pair(u, v);
                edge_t e;
                do {
                    e = edge_stack.top();
                    edge_stack.pop();
                } while(e != last_edge);
            }
        } else if(v != pred_vec[u] && d_vec[v] < d_vec[u]) {   // (u, v) is a back edge
            edge_stack.push(make_pair(u, v));
            low_vec[u] = min(low_vec[u], d_vec[v]);
        }
    }
}

//the single edge in a size 2 bcc is a bridge
void graph_t::find_bridge_edges(vector<edge_t>& out_vec)
{
    out_vec.clear();
    
    vector<graph_hash_t> bcc_vec;
    find_bicon_comp(bcc_vec);
    
    for(int i = 0; i < bcc_vec.size(); ++i) {
        if(bcc_vec[i].size() == 2) {
            edge_t e;
            e = *bcc_vec[i].edge_set.begin();
            out_vec.push_back(e);
        }
    }
}