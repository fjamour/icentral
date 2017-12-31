#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <stack>
#include <queue>
#include <numeric>
#include <cmath>

#include "bc.h"
#include "utility.h"

/*************************************************
 * EXPERIMENTAL
 *************************************************/
void Update_BC_mem(
        vector<double>&         BC_vec,
        graph_t&                graph,
        comp_type_t             comp_type,
        edge_t                  e,
        vector<iter_info_t>&    iter_info_vec
        )
{
    if(comp_type == GRAPH) {
        component_t comp;
        comp.comp_type = GRAPH;
        comp.subgraph.fill_graph(graph);
        vector<double> dBC_vec;
        iCentral_mem(dBC_vec, comp, e, iter_info_vec);
        for(int i = 0; i < dBC_vec.size(); ++i) {
            BC_vec[i] += dBC_vec[i];
        }
    }
}
void iCentral_mem(
        vector<double>&         dBC_vec,
        component_t&            comp,
        edge_t                  e,
        vector<iter_info_t>&    iter_info_vec
        )
{
    fill_vec<double>(dBC_vec, comp.subgraph.size(), 0.0);
    vector<int> d_src_vec, d_dst_vec;
    comp.subgraph.find_sssp(e.first, d_src_vec);
    comp.subgraph.find_sssp(e.second, d_dst_vec);    
    for(node_id_t s = 0; s < comp.subgraph.size(); ++s) {
        if(d_src_vec[s] != d_dst_vec[s]) {
            iCentral_iter_mem(dBC_vec, comp, s, e, iter_info_vec[s]);
        }

    }
}




//enum state_t {UP, DOWN, NT};
//void iCentral_iter_mem(
//        vector<double>& dBC_vec,    // delta BC of vertices
//        component_t&    comp,       // component could be BCC, MUC, or just a graph
//        node_id_t       s,          // source of the iteration
//        edge_t          e,          // inserted edge
//        iter_info_t&    iter_info   //TODO to be used later?
//        )
//{
//    subgraph_t& g = comp.subgraph;
//    
//    vector<int>&                 sigma_vec = iter_info.sigma_vec;
//    vector<int>&                 sigma_inc_vec = iter_info.sigma_inc_vec;
//    vector<int>&                 new_sigma_vec = iter_info.new_sigma_vec;
//    vector<int>&                 dist_map = iter_info.dist_vec;
//    vector<double>&              delta_vec = iter_info.delta_vec;
//    vector<double>&              new_delta_vec = iter_info.new_delta_vec;
//    
//    static vector<char>                 t_vec;
//    static vector<bool>                 moved_vec;
//    static vector<int>                  movement_vec;
//    static vector<queue<node_id_t> >    level_vec;
//    
//    queue<node_id_t>             Q;
//    
//    node_id_t src, dst;
//    if(dist_map[e.first] > dist_map[e.second]) {
//        src = e.second;
//        dst = e.first;
//    } else {
//        src = e.first;
//        dst = e.second;
//    }
//    
//    level_vec.resize(g.size() + 1);
//    fill_vec<char>(t_vec, g.size(), NT);
//    fill(sigma_inc_vec.begin(), sigma_inc_vec.end(), 0);
//    fill_vec<bool>(moved_vec, g.size(), false);
//    fill_vec<int>(movement_vec, g.size(), 0);
//    
//    sigma_inc_vec[dst] += sigma_vec[src];
//    t_vec[dst] = DOWN;
//    new_sigma_vec = sigma_vec;
//    if(dist_map[dst] != dist_map[src]+1)
//        new_sigma_vec[dst] = new_sigma_vec[src];
//    else
//        new_sigma_vec[dst] += new_sigma_vec[src];
//    moved_vec[dst] = true;
//    movement_vec[dst] = dist_map[dst] - dist_map[src];
//    level_vec[dist_map[dst]].push(dst);
//    Q.push(dst);
//
//    while(!Q.empty()) {
//        node_id_t v = Q.front(); Q.pop();
//        vector<node_id_t> nbr_vec = g.get_nbrs(v);
//        for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
//            node_id_t w = nbr_vec[nbr_idx];
//            int dist = (movement_vec[w]-movement_vec[v]) - (dist_map[v]-dist_map[w]+1);
//            if(t_vec[w] == NT && dist >= 0) {
//                t_vec[w] = DOWN;
//                Q.push(w);
//                if(dist > 0) {
//                    sigma_inc_vec[w] = sigma_inc_vec[v];
//                    new_sigma_vec[w] = sigma_inc_vec[v];
//                    moved_vec[w] = true;
//                    movement_vec[w] = dist;
//                } else if(dist == 0) {
//                    sigma_inc_vec[w] += sigma_inc_vec[v];
//                    new_sigma_vec[w] += sigma_inc_vec[v];
//                }
//            } else if(t_vec[w] == DOWN && dist >= 0) {
//                new_sigma_vec[w] += sigma_inc_vec[v];
//                sigma_inc_vec[w] += sigma_inc_vec[v];
//            }
//        }
//        dist_map[v] -= movement_vec[v];
//        level_vec[dist_map[v]].push(v);
//    }
//    
//    fill(new_delta_vec.begin(), new_delta_vec.end(), 0);
//    
//    for(int level = g.size(); level >= 0; --level) {
//        while(!level_vec[level].empty()) {
//            node_id_t w = level_vec[level].front(); level_vec[level].pop();
//            vector<node_id_t> nbr_vec = g.get_nbrs(w);
//            for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
//                node_id_t v = nbr_vec[nbr_idx];
//                if(dist_map[v] == dist_map[w]-1) {
//                  if(t_vec[v] == NT) {
//                      level_vec[level-1].push(v);
//                      t_vec[v] = UP;
//                      new_delta_vec[v] = delta_vec[v];
//                  }
//                  new_delta_vec[v] += ((double)new_sigma_vec[v]/new_sigma_vec[w])*(1+new_delta_vec[w]);
//                  if(t_vec[v] == UP && (v != src || w != dst)) {
//                      new_delta_vec[v] -= ((double)sigma_vec[v]/sigma_vec[w])*(1+delta_vec[w]);
//                  }
//                } else if(dist_map[v] == dist_map[w] && moved_vec[w] && !moved_vec[v]) {
//                    if(t_vec[v] == NT) {
//                        level_vec[level].push(v);
//                        t_vec[v] = UP;
//                        new_delta_vec[v] = delta_vec[v];
//                    }
//                }
//                if (w != s) {
//                    dBC_vec[w] += (new_delta_vec[w] - delta_vec[w]);
//                }
//            }
//        }
//    }
//    sigma_vec = new_sigma_vec;
//    for(node_id_t v = 0; v < g.size(); ++v) {
//        if(t_vec[v] != NT) {
//            delta_vec[v] = new_delta_vec[v];
//        }
//    }
//}



void iCentral_iter_mem(
        vector<double>& dBC_vec,    // delta BC of vertices
        component_t&    comp,       // component could be BCC, MUC, or just a graph
        node_id_t       s,          // source of the iteration
        edge_t          e,          // inserted edge
        iter_info_t&    iter_info   //TODO to be used later?
        )
{
//    RBFS(dBC_vec, comp, s, iter_info, false, true);
//    partial_BBFS(iter_info, comp, s, e);
//    RBFS(dBC_vec, comp, s, iter_info, true, false);
//    fill(iter_info.sigma_inc_vec.begin(), iter_info.sigma_inc_vec.end(), 0);
//    fill(iter_info.visited_vec.begin(), iter_info.visited_vec.end(), false);
    
    subgraph_t& g = comp.subgraph;
    
    vector<vector<node_id_t> >&  P = iter_info.P;
    vector<int>&                 sigma_vec = iter_info.sigma_vec;
    vector<int>&                 sigma_inc_vec = iter_info.sigma_inc_vec;
    vector<int>&                 dist_map = iter_info.dist_vec;
    vector<bool>&                visited_vec = iter_info.visited_vec;
    //vector<node_id_t>&           S = iter_info.S;
    
    vector<double>&              delta_vec = iter_info.delta_vec;
    
    vector<double>&              new_delta_vec = iter_info.new_delta_vec;
    new_delta_vec = delta_vec;
    
    static vector<vector<node_id_t> >   new_P;
    fill_vec<vector<node_id_t> >(new_P, g.size(), vector<node_id_t>());
    //static vector<vector<node_id_t> >   lost_P;;
    //fill_vec<vector<node_id_t> >(lost_P, g.size(), vector<node_id_t>());
    static vector<int>                  new_sigma_vec = iter_info.new_sigma_vec;
    new_sigma_vec = sigma_vec;
    static vector<stack<node_id_t> >    level_vec;

    queue<node_id_t>             Q;
    node_id_t src, dst;
    if(dist_map[e.first] > dist_map[e.second]) {
        src = e.second;
        dst = e.first;
    } else {
        src = e.first;
        dst = e.second;
    }
    
    //S.clear();
    level_vec.resize(dist_map[dst]+1);
    
    //Compute new path counts and paths
    if(dist_map[dst] != dist_map[src]+1) {
        //lost_P[dst] = P[dst];
        new_P[dst].push_back(src);
        for(int p = 0; p < P[dst].size(); ++p){
            node_id_t parent = P[dst][p];
            new_delta_vec[parent] -= (((double) sigma_vec[parent] / sigma_vec[dst]) * (1 + delta_vec[dst]));
            if(dist_map[parent] == dist_map[src]+1) {
                level_vec[dist_map[parent]].push(parent);
            }
        }
        P[dst].clear();
        //P[dst].push_back(src);
        new_sigma_vec[dst] = new_sigma_vec[src];
    } else {
        new_P[dst].push_back(src);
        //P[dst].push_back(src);
        new_sigma_vec[dst] += new_sigma_vec[src];
    }
        
    Q.push(dst);
    dist_map[dst] = dist_map[src] + 1;
    sigma_inc_vec[dst] = sigma_vec[src];
    visited_vec[dst] = true;

    while (!Q.empty()) {
        node_id_t v = Q.front();
        if(dist_map[v] >= level_vec.size()) {
            stack<node_id_t> qq;
            level_vec.push_back(qq);
            level_vec[dist_map[v]].push(v);
        } else {
            level_vec[dist_map[v]].push(v);
        }
        Q.pop();
        vector<node_id_t> nbr_vec = g.get_nbrs(v);
        for (int nbr_idx = 0; nbr_idx < nbr_vec.size(); ++nbr_idx) {
            node_id_t w = nbr_vec[nbr_idx];
            if (dist_map[w] > (dist_map[v] + 1)) {
                dist_map[w] = dist_map[v] + 1;
                //lost_P[w] = P[w];
                new_P[w].push_back(v);
                for(int p = 0; p < P[w].size(); ++p){
                    node_id_t parent = P[w][p];
                    new_delta_vec[parent] -= (((double) sigma_vec[parent] / sigma_vec[w]) * (1 + delta_vec[w]));
                }
                P[w].clear();
                //P[w].push_back(v);
                new_sigma_vec[w] = 0;
                sigma_inc_vec[w] = sigma_inc_vec[v];
                new_sigma_vec[w] += sigma_inc_vec[w];
                if (!visited_vec[w]) {
                    visited_vec[w] = true;
                    //if(new_sigma_vec[w] != sigma_vec[w])
                        Q.push(w);
                }
            } else if (dist_map[w] == (dist_map[v] + 1)) {
                sigma_inc_vec[w] += sigma_inc_vec[v];
                new_sigma_vec[w] += sigma_inc_vec[v];
                if (find(P[w].begin(), P[w].end(), v) == P[w].end()) {
                    new_P[w].push_back(v);
                }
                if (!visited_vec[w]) {
                    visited_vec[w] = true;
                    //if(new_sigma_vec[w] != sigma_vec[w])
                        Q.push(w);
                }
            }
        }
    }
    
    fill_vec<bool>(visited_vec, g.size(), false);
    
    for(int level = level_vec.size()-1; level >= 0; --level) {
        while(!level_vec[level].empty()) {
            node_id_t w = level_vec[level].top();
            level_vec[level].pop();
            if(!visited_vec[w]) {
                visited_vec[w] = true;
                for (int idx = 0; idx < P[w].size(); ++idx) {
                    node_id_t v = P[w][idx];
                    new_delta_vec[v] -= (((double) sigma_vec[v] / sigma_vec[w]) * (1 + delta_vec[w]));
                    new_delta_vec[v] += (((double) new_sigma_vec[v] / new_sigma_vec[w]) * (1 + new_delta_vec[w]));
                    if(level > 0)
                        level_vec[level-1].push(v);
                }
                for (int idx = 0; idx < new_P[w].size(); ++idx) {
                    node_id_t v = new_P[w][idx];
                    new_delta_vec[v] += (((double) new_sigma_vec[v] / new_sigma_vec[w]) * (1 + new_delta_vec[w]));
                    if(level > 0)
                        level_vec[level-1].push(v);
                }
                //has to update P with new_P for the next edge insertions
                P[w] = new_P[w];
                if (w != s) {
                    dBC_vec[w] -= delta_vec[w]/2.0;
                    dBC_vec[w] += new_delta_vec[w]/2.0;
                }
            }
        }
    }
    
    iter_info.delta_vec = new_delta_vec;
    iter_info.sigma_vec = new_sigma_vec;
    fill(iter_info.sigma_inc_vec.begin(), iter_info.sigma_inc_vec.end(), 0);
    fill(iter_info.visited_vec.begin(), iter_info.visited_vec.end(), false);
}
