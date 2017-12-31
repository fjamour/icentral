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
#include "mcb_find.h"
#include "bc.h"
#include "utility.h"

using namespace std;


vector<node_id_t>& graph_t::get_nbrs(node_id_t id)
{
    return nodes_vec[id].nbrs_vec;
}
    
void graph_t::insert_edge(node_id_t src, node_id_t dst)
{
    pair<node_id_t, node_id_t> edge = make_pair(src, dst);
    if(edge_set.find(edge) != edge_set.end())
        return;
    edge = make_pair(dst, src);
    if(edge_set.find(edge) != edge_set.end())
        return;
    edge_set.insert(edge);
    nodes_vec[src].nbrs_vec.push_back(dst);
    nodes_vec[dst].nbrs_vec.push_back(src);
}

void graph_t::remove_edge(node_id_t src, node_id_t dst)
{
    pair<node_id_t, node_id_t> edge = make_pair(src, dst);
    edge_set.erase(edge);
    edge = make_pair(dst, src);
    edge_set.erase(edge);
    
    vector<node_id_t>& v1 = nodes_vec[src].nbrs_vec;
    vector<node_id_t>& v2 = nodes_vec[dst].nbrs_vec;
    v1.erase(std::remove(v1.begin(), v1.end(), dst), v1.end());
    v2.erase(std::remove(v2.begin(), v2.end(), src), v2.end());
}

void graph_t::init_size(size_t num_nodes)
{
    nodes_vec.resize(num_nodes);
    node_to_muc_vec.resize(num_nodes);
    fill(node_to_muc_vec.begin(), node_to_muc_vec.end(), INF);
    i_init_internals();
}

void graph_t::i_init_internals()
{
    tmp_bool_vec.resize(size());
    //bc_vec.resize(size());
    bc_computed = false;
}

size_t graph_t::size()
{
    return nodes_vec.size();
}

/*
 * #vertices    #edges
 * src1         dst1
 * .
 * .
 * src#edges    dst#edges
 */
void graph_t::read_graph(string path)
{
    ifstream fin;
    fin.open(path.c_str(), ios::in);
    if(!fin.good()) {
        printf("Can't open the file [%s]\n", path.c_str());
        exit(1);
    }
    size_t N, M;
    fin >> N >> M;
    init_size(N);
    for(int i = 0; i < M; ++i) {
        node_id_t src, dst;
        fin >> src >> dst;
        //TODO could this cause issues?
        if(src < N && dst < N) {
            insert_edge(src, dst);
        }
    }
    fin.close();
    i_init_internals();
}

void graph_t::print_edgelist()
{
     for(set<pair<node_id_t, node_id_t> >::iterator it = edge_set.begin(); it != edge_set.end(); it++) {
        printf("[%d] === [%d]\n", it->first, it->second);
    }   
}

void graph_t::print_mcbs()
{
    for(int c = 0; c < mcb.cycle_vec.size(); ++c) {
        printf("Cycle [%d]:\n", c);
        for(int e = 0; e < mcb.cycle_vec[c].size(); ++e) {
            printf("    [%d]===[%d]\n", mcb.cycle_vec[c][e].first, mcb.cycle_vec[c][e].second);
        }        
    }
}

void graph_t::print_node_to_muc()
{
    for(int i = 0; i < node_to_muc_vec.size(); ++i) {
       printf("[%d] --> MUC [%d]\n", i, node_to_muc_vec[i]); 
    }
}

void graph_t::get_shortest_path(
        node_id_t src,
        node_id_t dst,
        vector<node_id_t>& out_path)
{
    //do a BFS from dst to src and store path
    vector<int>         dist_vec;
    vector<node_id_t>   parent_vec;
    dist_vec.resize(size());
    parent_vec.resize(size());
    fill(dist_vec.begin(), dist_vec.end(), -1);
    fill(parent_vec.begin(), parent_vec.end(), -1);
    
    queue<node_id_t> q;
    dist_vec[dst] = 0;
    q.push(dst);
    while(!q.empty()) {
        node_id_t node = q.front(); q.pop();
        if(node == src)
            break;
        vector<node_id_t> nbr_vec = get_nbrs(node);
        for(int i = 0; i < nbr_vec.size(); ++i) {
            node_id_t nbr_id = nbr_vec[i];
            if(dist_vec[nbr_id] == -1) {
                dist_vec[nbr_id] = dist_vec[node] + 1;
                parent_vec[nbr_id] = node;
                q.push(nbr_id);
            }
        }
    }
    //fill the output path vector:
    out_path.clear();
    node_id_t nd = src;
    while(parent_vec[nd] != -1) {
        out_path.push_back(nd);
        nd = parent_vec[nd];
    }
    out_path.push_back(nd);//insert dst
}

bool graph_t::has_edge(const edge_t& e)
{
    edge_t e1 = e;
    edge_t e2; e2.first = e1.second; e2.second = e1.first;
    return (edge_set.find(e1) != edge_set.end() ||
            edge_set.find(e2) != edge_set.end());
}

//assume @out_vec has same size as graph and all elements initialized to -1
void graph_t::find_sssp(node_id_t s, vector<int>& out_vec)
{
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

void graph_t::insert_edge_update_bc_experimental(node_id_t src,
                          node_id_t dst,
                          bcc_stat_t& bcc_stat,
                          int max_iter_d1,
                          int max_iter_d2)
{
    /*
     * 1. edge must be in the graph so that bcc extraction works properly
     * 2. the subgraph in the bcc must not have the edge (done inside bcc_delta_t::compute_bc)
     */
    timer tm;//
    
    insert_edge(src, dst);    
    bcc_delta_t bcc;
    tm.start();//
    find_edge_bcc(bcc, src, dst);
    tm.stop();//
    bcc_stat.bcc_num_nodes = bcc.bcc_subgraph.nodes_map.size();
    bcc_stat.bcc_num_edges = bcc.bcc_subgraph.edge_set.size();
    bcc_stat.bcc_find_time = tm.interval();//
    tm.start();
    if(max_iter_d1 != -1 || max_iter_d2 != -1)
        bcc.compute_bc_maxiter_exp(bc_map, src, dst, bcc_stat, max_iter_d1, max_iter_d2);
    else
        bcc.compute_bc_exp(bc_map, src, dst, bcc_stat, max_iter_d1, max_iter_d2);
    tm.stop();//
    bcc_stat.bc_update_time = tm.interval();//
    
    bcc.bcc_subgraph.remove_edge(src, dst);
    bcc.bcc_subgraph.find_pruning_counts_exp(src, 
                        dst, bcc_stat.tot_d0_iter, 
                        bcc_stat.tot_d1_iter, 
                        bcc_stat.tot_d2_iter);
    bcc.bcc_subgraph.insert_edge(src, dst);
}

void graph_t::approx_bcc_iter_tm(node_id_t src,
                                 node_id_t dst,
                                 double& avg_iter_time,
                                 int num_iter)
{
    timer tm;
    
    bcc_scratch_t bcc;
    find_edge_bcc(bcc, src, dst);
    if(num_iter == -1) {
        num_iter = bcc.bcc_subgraph.size();
    }
    
    tr1_map_t(double) bc_map;
    bcc.bcc_subgraph.i_fill_map<double>(bc_map, 0);
    tm.start();
    bcc.compute_bc(bc_map, num_iter);
    tm.stop();
    
    avg_iter_time = tm.interval()/num_iter;
}

void graph_t::init_bc()
{
    if(!bc_computed) {
        bc_vec = brandes_bc(*this);
        for(int i = 0; i < bc_vec.size(); ++i) {
            bc_map.insert(make_pair(i, bc_vec[i]));
        }
        bc_computed = true;
    }    
}

void graph_t::tmp_fun()
{
    print_edgelist();
    printf("===================\n");
    print_mucs(false);
}

/************************************************************/
/************************************************************/

//IMP: assumes @e is not in the graph already
void graph_t::find_edge_bcc(
                component_t&    comp,
                edge_t          e
                )
{
    //a. find the bcc subgraph
    //b. find the articulation points in bcc
    //c. find the sizes of the subgraphs connected to each art point
    if(has_edge(e)) {
        printf("WARNING: inserted duplicate edge!\n");
    }
    insert_edge(e.first, e.second);
    graph_hash_t g;
    find_edge_bcc_subgraph(g, e.first, e.second);
    //g.remove_edge(e.first, e.second); //TODO: make sure it's okay! [NO, needed below]
    vector<node_id_t> art_pt_vec;
    find_art_points(art_pt_vec);
    set<node_id_t> art_pt_set;//to make art points unique
    for(int i = 0; i < art_pt_vec.size(); ++i) {
        art_pt_set.insert(art_pt_vec[i]);
    }
    
    graph_hash_t::nodes_map_t::iterator it;
    for(it = g.nodes_map.begin();
        it != g.nodes_map.end();
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
                if(!g.has_edge(v, u) && !visited_vec[u]) {
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
                    //printf("---%d, %d", v, cnt);
                }
            }
            comp.art_pt_map.insert(make_pair(v, subgraph_sz_vec));
        }
    }
    //fix art_pt_map
    g.remove_edge(e.first, e.second);
    comp.subgraph.fill_graph(g);
    
    component_t::art_pt_map_t new_art_pt_map;
    for(component_t::art_pt_map_t::iterator it = comp.art_pt_map.begin();
            it != comp.art_pt_map.end();
            ++it) {
        node_id_t new_v = comp.subgraph.outin_label_map[it->first];
        new_art_pt_map.insert(make_pair(new_v, it->second));
    }
    comp.art_pt_map = new_art_pt_map;
    
    remove_edge(e.first, e.second);
}