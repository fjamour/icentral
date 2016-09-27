#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include <stack>

#include <numeric>

#include "graph_t.h"
#include "mcb_find.h"
#include "utility.h"

muc_t::muc_t()
{
    valid = true;
}

//returns true is muc and cycle share a vertex
bool shared_vertex(muc_t& muc, cycle_t& cycle)
{
    bool res = false;
    for(int i = 0; i < cycle.size(); ++i) {
        if(   muc.muc_subgraph.nodes_map.find(cycle[i].first) != muc.muc_subgraph.nodes_map.end() 
           || muc.muc_subgraph.nodes_map.find(cycle[i].second) != muc.muc_subgraph.nodes_map.end())
            return true;
    }
    return res;
}

bool shared_vertex(cycle_t& c1, cycle_t& c2)
{
   bool res = false;
    for(int i = 0; i < c1.size(); ++i) {
        for(int j = 0; j < c2.size(); ++j) {
            if(   c1[i].first   == c2[j].first
               || c1[i].first   == c2[j].second
               || c1[i].second  == c2[j].first
               || c1[i].second  == c2[j].second) {
                return true;
            }
        }
    }
    return res;    
}

void graph_t::find_mucs_fast()
{
    //1. find bridges
    //2. delete bridges
    //3. find connected components of size >=3 (these are MUCs)
    vector<edge_t> bridge_vec;
    find_bridge_edges(bridge_vec);
    graph_hash_t gh;
//    printf("DBG filling graph hash...\n");
    gh.fill_graph(*this);
//    printf("DBG removing bridges...\n");
    for(int i = 0; i < bridge_vec.size(); ++i) {
        edge_t e = bridge_vec[i];
        gh.remove_edge(e.first, e.second);
    }
    vector<graph_hash_t> conn_comp_vec;
//    printf("DBG finding connected components...\n");
    gh.find_conn_comp(conn_comp_vec);
//    printf("DBG constructing MUCs...\n");
    for(int i = 0; i < conn_comp_vec.size(); ++i) {
        if(conn_comp_vec[i].size() >= 3) {
            muc_t muc;
            muc.id = muc_vec.size();
            muc.muc_subgraph = conn_comp_vec[i];
            graph_hash_t::nodes_map_t::iterator it;
            it = muc.muc_subgraph.nodes_map.begin();
            for(; it != muc.muc_subgraph.nodes_map.end(); ++it) {
                node_id_t v = it->first;
                node_to_muc_vec[v] = muc.id;
            }
            muc_vec.push_back(muc);
        }
    }
//    printf("DBG constructing single node MUCs...\n");
    //print_mucs();
    //create single vertex muc's
    for(node_id_t i = 0; i < node_to_muc_vec.size(); ++i) {
        if(node_to_muc_vec[i] == INF) {
            muc_t muc;
            muc.id = muc_vec.size();
            muc.muc_subgraph.insert_node(i);
            node_to_muc_vec[i] = muc.id;
            muc_vec.push_back(muc);
        }
    }
//    printf("DBG finding connection vertices...\n");
    find_conn_verts();
    
    //find the disconnected subgraphs that will result from the deletion of
    //connection vertices
//    printf("DBG finding MUC subgraphs...\n");
    for(muc_id_t m = 0; m < muc_vec.size(); ++m) {
        find_muc_subgraphs(m);
    }
//    printf("DBG DONE!...\n");
}

void graph_t::find_muc_mcb()
{
    mcb_find(this, &mcb);
    vector<bool> visited_vec;
    visited_vec.resize(mcb.cycle_vec.size());
    fill(visited_vec.begin(), visited_vec.end(), false);
    graph_hash_t cycle_graph;
    for(int i = 0; i < mcb.cycle_vec.size(); ++i) {
        cycle_graph.insert_node(i);
    }
    for(int i = 0; i < mcb.cycle_vec.size(); ++i) {
        for(int j = 0; j < mcb.cycle_vec.size(); ++j) {
            if(i != j && shared_vertex(mcb.cycle_vec[i], mcb.cycle_vec[j])) {
                cycle_graph.insert_edge(i, j);
            }
        }        
    }
    //cycle_graph.print_graph();
    for(int i = 0; i < mcb.cycle_vec.size(); ++i) {
        if(!visited_vec[i]) {
            muc_t muc;
            muc.id = muc_vec.size();
            merge_muc_cycle(muc, mcb.cycle_vec[i]);       
            int src = i;
            visited_vec[src] = true;
            queue<int> q;
            q.push(src);
            while(!q.empty()) {
                int node = q.front(); q.pop();
                for(int nb_i = 0; nb_i < cycle_graph.nodes_map[node].nbrs_vec.size(); ++nb_i) {
                    int nbr = cycle_graph.nodes_map[node].nbrs_vec[nb_i];
                    if(!visited_vec[nbr]) {
                        visited_vec[nbr] = true;
                        merge_muc_cycle(muc, mcb.cycle_vec[nbr]);
                        q.push(nbr);
                    }
                }
            }
            muc_vec.push_back(muc);
        }
    }
    //print_mucs();
    //create single vertex muc's
    for(node_id_t i = 0; i < node_to_muc_vec.size(); ++i) {
        if(node_to_muc_vec[i] == INF) {
            muc_t muc;
            muc.id = muc_vec.size();
            muc.muc_subgraph.insert_node(i);
            node_to_muc_vec[i] = muc.id;
            muc_vec.push_back(muc);
        }
    }
    
    find_conn_verts();
    
    //find the disconnected subgraphs that will result from the deletion of
    //connection vertices
    for(muc_id_t m = 0; m < muc_vec.size(); ++m) {
        find_muc_subgraphs(m);
    }
}

//add all the edges in cycle to muc
void graph_t::merge_muc_cycle(muc_t& muc, cycle_t& cycle)
{
    for(int i = 0; i < cycle.size(); ++i) {
        node_id_t src, dst;
        src = cycle[i].first;
        dst = cycle[i].second;
        muc.muc_subgraph.insert_edge(src, dst);
        node_to_muc_vec[src] = muc.id;
        node_to_muc_vec[dst] = muc.id;
    }
}

void graph_t::find_mucs()
{   
    //find_muc_mcb();
    find_mucs_fast();
}

void graph_t::find_conn_verts()
{
    for(muc_id_t m = 0; m < muc_vec.size(); ++m) {
        if(muc_vec[m].valid) {
            muc_vec[m].conn_vertex_map.clear();
            muc_vec[m].subgraph_map.clear();
        }
    }
    //find the connection vertices for each muc
    //edges that belong to the graph and don't belong
    //to any muc connect two connection vertices
    //1. collect all edges in mucs in one set
    //2. do simple subtraction of edges and update muc's
    set<edge_t> all_muc_edge_set;
    for(muc_id_t i = 0; i < muc_vec.size(); ++i) {
        if(muc_vec[i].valid) {
            set<edge_t>::iterator it1 = muc_vec[i].muc_subgraph.edge_set.begin();
            set<edge_t>::iterator it2 = muc_vec[i].muc_subgraph.edge_set.end();
            all_muc_edge_set.insert(it1, it2);        
        }
    }
    for(set<edge_t>::iterator it = edge_set.begin(); it != edge_set.end(); ++it) {
        edge_t edge   = *it;
        edge_t edge_r = make_pair(edge.second, edge.first);
        if(all_muc_edge_set.find(edge) == all_muc_edge_set.end() && all_muc_edge_set.find(edge_r) == all_muc_edge_set.end()) {//edge not in muc's, connects two connection vertices
            edge_t bridge_edge = *it;
            muc_id_t muc1_id = node_to_muc_vec[bridge_edge.first];
            muc_id_t muc2_id = node_to_muc_vec[bridge_edge.second];
            muc_vec[muc1_id].insert_conn_vertex(bridge_edge.first, bridge_edge.second);
            muc_vec[muc2_id].insert_conn_vertex(bridge_edge.second, bridge_edge.first);
        }
    }   
}

//TODO  not sure what to do with single node muc's, will not
//      do anything for them now
void graph_t::find_muc_subgraphs(muc_id_t muc_id) {
    //THIS IS SLOW, and not needed for the QUBE ideal speedup experiemtn
    //so just return if you want to do this experiment..
    //return;
    if (muc_vec[muc_id].muc_subgraph.size() == 1)
        return;
    if(!muc_vec[muc_id].valid)
        return;
    for (muc_t::conn_vert_map_t::iterator
        it = muc_vec[muc_id].conn_vertex_map.begin();
            it != muc_vec[muc_id].conn_vertex_map.end();
            ++it) {
        node_id_t conn_vert = it->first;
        graph_hash_t g;
        for (vector<node_id_t>::iterator it2 = it->second.begin();
                it2 != it->second.end(); ++it2) {
            node_id_t bfs_source = *it2;
            vector<bool> visited_vec;
            visited_vec.resize(size());
            fill(visited_vec.begin(), visited_vec.end(), false);
            visited_vec[bfs_source] = true;
            visited_vec[conn_vert] = true;
            queue<node_id_t> q;
            q.push(bfs_source);
            while (!q.empty()) {
                node_id_t node = q.front();
                q.pop();
                g.insert_node(node);
                for (int bfs_nbr_i = 0; bfs_nbr_i < get_nbrs(node).size(); ++bfs_nbr_i) {
                    node_id_t bfs_nbr = get_nbrs(node)[bfs_nbr_i];
                    if (!visited_vec[bfs_nbr]) {
                        q.push(bfs_nbr);
                        visited_vec[bfs_nbr] = true;
                        if (bfs_nbr != conn_vert) {
                            g.insert_edge(node, bfs_nbr);
                        }
                    }
                }
            }
        }
        muc_vec[muc_id].subgraph_map.insert(make_pair(conn_vert, g));
    }
}

void graph_t::print_mucs(bool with_nodes)
{
    for(int i = 0; i < muc_vec.size(); ++i) {
        if(!muc_vec[i].valid)
            continue;
        printf("MUC [%d] ====================================\n", i);
        printf("Connection vertices: ");
        for(muc_t::conn_vert_map_t::iterator
            it =  muc_vec[i].conn_vertex_map.begin();
            it != muc_vec[i].conn_vertex_map.end();
            ++it) {
            printf("[%d |G|=%d] ", it->first, muc_vec[i].subgraph_map[it->first].size());
        }
        printf("\n");
        printf("Bridges: ");
        for(muc_t::conn_vert_map_t::iterator
            it =  muc_vec[i].conn_vertex_map.begin();
            it != muc_vec[i].conn_vertex_map.end();
            ++it) {
            for(vector<node_id_t>::iterator it2 = it->second.begin();
                it2 != it->second.end(); ++it2) {
                printf("[%d == %d] ", it->first, *it2);
            }
            
        }
        printf("\n");
        muc_vec[i].muc_subgraph.print_graph(with_nodes);
    }
}

//make sure the node with @id exists in the muc!
//this function will blindly insert it
void muc_t::insert_conn_vertex(node_id_t id, node_id_t nbr) {
    if(conn_vertex_map.find(id) != conn_vertex_map.end()) {
        conn_vertex_map[id].push_back(nbr);
    } else {
        vector<node_id_t> vec;
        vec.push_back(nbr);
        conn_vertex_map.insert(make_pair(id, vec));
    }
}

//muc1 = muc1 U muc2
void graph_t::merge_muc_muc(muc_t& muc1, muc_t& muc2)
{
    //1. copy the graph of muc2 into muc1
    //2. for each bridge edge, if the other end becomes in the graph, add the edge
    //3. mark muc2 as invalid
    //4. ASSUME muc1.id is valid, and update node_to_muc_vec for vertices in muc2
    //   to hold the id of muc1
    muc1.muc_subgraph.nodes_map.insert(
                muc2.muc_subgraph.nodes_map.begin(),
                muc2.muc_subgraph.nodes_map.end());
    muc1.muc_subgraph.edge_set.insert(
                muc2.muc_subgraph.edge_set.begin(),
                muc2.muc_subgraph.edge_set.end());
    for(muc_t::conn_vert_map_t::iterator it = muc2.conn_vertex_map.begin();
                                         it != muc2.conn_vertex_map.end(); 
                                         ++it) {
        node_id_t conn_vert = it->first;
        for(int i = 0; i < it->second.size(); ++i) {
            if(muc1.muc_subgraph.nodes_map.find(it->second[i]) != muc1.muc_subgraph.nodes_map.end()) {
                muc1.muc_subgraph.insert_edge(conn_vert, it->second[i]);
            }
        }
    }
    muc2.valid = false;
    for(graph_hash_t::nodes_map_t::iterator it = muc2.muc_subgraph.nodes_map.begin();
                                            it != muc2.muc_subgraph.nodes_map.end();
                                            ++it) {
        node_to_muc_vec[it->first] = muc1.id;
    }
}


void graph_t::insert_edge_update_muc(node_id_t src, node_id_t dst)
{
    muc_id_t src_muc_id = node_to_muc_vec[src];
    muc_id_t dst_muc_id = node_to_muc_vec[dst];
    if(src_muc_id == dst_muc_id) {
        //case 1: both @src and @dst are in the same muc
        //        then just update the graph and the muc graph
        //just update the muc and the graph
        muc_vec[src_muc_id].muc_subgraph.insert_edge(src, dst);
        insert_edge(src, dst);
    } else {
        //case 2: @src and @dst are not in the same muc
        muc_t new_muc;
        vector<node_id_t> shortest_path;
        get_shortest_path(src, dst, shortest_path);
        new_muc.id = muc_vec.size();
        merge_muc_muc(new_muc, muc_vec[dst_muc_id]);
        merge_muc_muc(new_muc, muc_vec[src_muc_id]);
        for(int i = 0; i < shortest_path.size(); ++i) {
            muc_id_t node_muc_id = node_to_muc_vec[shortest_path[i]];
            if(node_muc_id != new_muc.id) {
                merge_muc_muc(new_muc, muc_vec[node_muc_id]);
            }
        }
        new_muc.muc_subgraph.insert_edge(src, dst);
        //cout << endl << endl;
        //new_muc.muc_subgraph.print_graph();
        muc_vec.push_back(new_muc);
        insert_edge(src, dst);
        //TODO updating the connection vertices must be done in a much faster way!
        find_conn_verts();
        find_muc_subgraphs(new_muc.id);
    }   
}

size_t graph_t::get_num_mucs()
{
    size_t res = 0;
    for(int i = 0; i < muc_vec.size(); ++i) {
        if(muc_vec[i].valid) {
            ++res;
        }
    }
    return res;
}




//void muc_t::compute_bc(tr1::unordered_map<node_id_t, double>& bc_map)
//{
//
//    //bc_map will have for each vertex in the muc it's betweenness centrality
//    //init bc of all nodes to zero
//    bc_map.clear();
//    if(!valid)
//        return;
//    for(graph_hash_t::nodes_map_t::iterator
//                        it =  muc_subgraph.nodes_map.begin();
//                        it != muc_subgraph.nodes_map.end();
//                        ++it) {
//        bc_map.insert(make_pair(it->first, 0));
//    }
//    
//
//    //do BFS's from the nodes in the muc
//    for(graph_hash_t::nodes_map_t::iterator
//                        it =  muc_subgraph.nodes_map.begin();
//                        it != muc_subgraph.nodes_map.end();
//                        ++it) {
//        node_id_t s = it->first;
//        
//        stack<node_id_t>                        S;
//        tr1_map_t(vector<node_id_t> )&          P = muc_subgraph.P; 
//        tr1_map_t(int)&                         path_cnt_map = muc_subgraph.path_cnt_vec;
//        tr1_map_t(int)&                         dist_map = muc_subgraph.dist_vec;
//        tr1_map_t(double)&                      pair_dep_map = muc_subgraph.pair_dep_vec;
//        tr1_map_t(double)&                      sigma_t_map = muc_subgraph.sigma_t_map; 
//        queue<node_id_t>&                       Q = muc_subgraph.Q;
//        
//        muc_subgraph.init_maps();           
//        
//        path_cnt_map[s] = 1;
//        dist_map[s] = 0;
//        
//        Q.push(s);
//        while(!Q.empty()) {
//            node_id_t v_i = Q.front(); Q.pop();
//            S.push(v_i);
//            for(int i = 0; i < muc_subgraph.nodes_map[v_i].nbrs_vec.size(); ++i) {
//                node_id_t v_n = muc_subgraph.nodes_map[v_i].nbrs_vec[i];
//                if(dist_map[v_n] < 0) {
//                    Q.push(v_n);
//                    dist_map[v_n] = dist_map[v_i] + 1;
//                }
//                if(dist_map[v_n] == dist_map[v_i] + 1) {
//                    path_cnt_map[v_n] = path_cnt_map[v_n] + path_cnt_map[v_i];
//                    P[v_n].push_back(v_i);
//                }
//            }
//        }
//        while(!S.empty()) {
//            node_id_t v_n = S.top(); S.pop();
//            if(   conn_vertex_map.find(s) != conn_vertex_map.end()
//               && conn_vertex_map.find(v_n) != conn_vertex_map.end()
//               && s != v_n) {
//                int VG_s, VG_n;
//                VG_s = subgraph_map[s].size();
//                VG_n = subgraph_map[v_n].size();
//                int c_t = VG_s*VG_n;
//                sigma_t_map[v_n] =  sigma_t_map[v_n] + c_t;
//                bc_map[v_n] = bc_map[v_n] + c_t;
//            }
//            for(int i = 0; i < P[v_n].size(); ++i) {
//                node_id_t v_p = P[v_n][i];
//                double sp_sn = ((double)path_cnt_map[v_p]/path_cnt_map[v_n]);
//                pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn*(1+pair_dep_map[v_n]);
//                if(conn_vertex_map.find(s) != conn_vertex_map.end()) {
//                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
//                    bc_map[v_p] = bc_map[v_p] + sigma_t_map[v_n]*sp_sn;
//                }
//            }
//            if(s != v_n) {
//                bc_map[v_n] = bc_map[v_n] + pair_dep_map[v_n];
//            }
//            if(conn_vertex_map.find(s) != conn_vertex_map.end()) {
//                int VG_s = subgraph_map[s].size();
//                bc_map[v_n] = bc_map[v_n] + pair_dep_map[v_n]*VG_s*2;
//            }
//        }
//    }
//    for(subgraph_map_t::iterator
//        it =  subgraph_map.begin();
//        it != subgraph_map.end();
//        ++it) {
//        vector<int> size_vec;
//        it->second.conn_comp_sizes(size_vec);
//        if(size_vec.size() > 1) {
//            int sub = 0;
//            for(int i = 0; i < size_vec.size(); ++i) {
//                sub += size_vec[i]*size_vec[i];
//            }
//            int VG_i = it->second.size();
//            bc_map[it->first] = bc_map[it->first] + VG_i*VG_i - sub;
//        }
//    }
//    
//    for(tr1::unordered_map<node_id_t, double>::iterator
//            it =  bc_map.begin();
//            it != bc_map.end();
//            ++it) {
//        it->second = it->second/2;
//    }
//}

void muc_t::compute_bc(tr1::unordered_map<node_id_t, double>& bc_map,
        int max_iter)
{

    //bc_map will have for each vertex in the muc it's betweenness centrality
    //init bc of all nodes to zero
    bc_map.clear();
    if(!valid)
        return;
    
    muc_fast_subgraph.fill_graph(muc_subgraph);
    tmp_conn_vertex_map.clear();
    tmp_subgraph_map.clear();
    for(conn_vert_map_t::iterator it = conn_vertex_map.begin();
            it != conn_vertex_map.end();
            ++it) {
        node_id_t new_v = muc_fast_subgraph.outin_label_map[it->first];
        tmp_conn_vertex_map.insert(make_pair(new_v, it->second));
    }
    for(subgraph_map_t::iterator it = subgraph_map.begin();
            it != subgraph_map.end();
            ++it) {
        node_id_t new_v = muc_fast_subgraph.outin_label_map[it->first];
        tmp_subgraph_map.insert(make_pair(new_v, it->second));
    }
 
    
    for(graph_hash_t::nodes_map_t::iterator
                        it =  muc_subgraph.nodes_map.begin();
                        it != muc_subgraph.nodes_map.end();
                        ++it) {
        bc_map.insert(make_pair(it->first, 0));
    }
    

    //do BFS's from the nodes in the muc
    if(max_iter == -1)
        max_iter = muc_fast_subgraph.nodes_vec.size();
    else
        max_iter = min(max_iter, (int)muc_fast_subgraph.nodes_vec.size());
    for(int i = 0; i < max_iter; ++i) {
        node_id_t s = i;
        
#ifdef USE_OPTIMIZED_BRANDES_QUBE
        stack<node_id_t>                    S;
        vector<vector<node_id_t> >&         P = muc_fast_subgraph.P; 
        vector<int>&                        path_cnt_map = muc_fast_subgraph.path_cnt_vec;
        vector<int>&                        dist_map = muc_fast_subgraph.dist_vec;
        vector<double>&                     pair_dep_map = muc_fast_subgraph.pair_dep_vec;
        vector<double>&                     sigma_t_map = muc_fast_subgraph.sigma_t_map; 
        queue<node_id_t>&                   Q = muc_fast_subgraph.Q;
        muc_fast_subgraph.init_maps();
#else
        stack<node_id_t>                   S;
        vector<vector<node_id_t> >         P;
        vector<int>                        path_cnt_map;
        vector<int>                        dist_map;
        vector<double>                     pair_dep_map;
        vector<double>                     sigma_t_map;
        queue<node_id_t>                   Q;
        int N = muc_fast_subgraph.size();
        fill_vec<vector<node_id_t> >(P, N, vector<node_id_t>());
        fill_vec<int>(path_cnt_map, N, 0);
        fill_vec<int>(dist_map, N, -1);
        fill_vec<double>(pair_dep_map, N, 0);
        fill_vec<double>(sigma_t_map, N, 0);
#endif
        
        path_cnt_map[s] = 1;
        dist_map[s] = 0;
        
        Q.push(s);
        while(!Q.empty()) {
            node_id_t v_i = Q.front(); Q.pop();
            S.push(v_i);
            for(int i = 0; i < muc_fast_subgraph.nodes_vec[v_i].size(); ++i) {
                node_id_t v_n = muc_fast_subgraph.nodes_vec[v_i][i];
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
        while(!S.empty()) {
            node_id_t v_n = S.top(); S.pop();
            if(   tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()
               && tmp_conn_vertex_map.find(v_n) != tmp_conn_vertex_map.end()
               && s != v_n) {
                int VG_s, VG_n;
                VG_s = tmp_subgraph_map[s].size();
                VG_n = tmp_subgraph_map[v_n].size();
                int c_t = VG_s*VG_n;
                sigma_t_map[v_n] =  sigma_t_map[v_n] + c_t;
                
                node_id_t new_v_n = muc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] + c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)path_cnt_map[v_p]/path_cnt_map[v_n]);
                pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn*(1+pair_dep_map[v_n]);
                if(tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()) {
                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
                    
                    node_id_t new_v_p = muc_fast_subgraph.inout_label_map[v_p];
                    bc_map[new_v_p] = bc_map[new_v_p] + sigma_t_map[v_n]*sp_sn;
                }
            }
            if(s != v_n) {
                node_id_t new_v_n = muc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] + pair_dep_map[v_n];
            }
            if(tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()) {
                int VG_s = tmp_subgraph_map[s].size();
                node_id_t new_v_n = muc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] + pair_dep_map[v_n]*VG_s*2;
            }
        }
    }
    for(subgraph_map_t::iterator
        it =  tmp_subgraph_map.begin();
        it != tmp_subgraph_map.end();
        ++it) {
        vector<int> size_vec;
        it->second.conn_comp_sizes(size_vec);
        if(size_vec.size() > 1) {
            int sub = 0;
            for(int i = 0; i < size_vec.size(); ++i) {
                sub += size_vec[i]*size_vec[i];
            }
            int VG_i = it->second.size();
            
            node_id_t new_node_id = muc_fast_subgraph.inout_label_map[it->first];
            bc_map[new_node_id] = bc_map[new_node_id] + VG_i*VG_i - sub;
        }
    }
    
    for(tr1::unordered_map<node_id_t, double>::iterator
            it =  bc_map.begin();
            it != bc_map.end();
            ++it) {
        it->second = it->second/2;
    }
}




/********************************************
* Incremental algorithm on top of MUC
* 
 *********************************************/

//////////////////////////////////////////////////
//IMP:: FOR NOW THIS GUY DOESN'T USE fast_subgraph
//////////////////////////////////////////////////
//void muc_t::compute_bc_inc(tr1_map_t(double)& bc_map,
//                        node_id_t src,
//                        node_id_t dst)
//{
//    if(!valid)
//        return;
//
//    //bc_map will have for each vertex in the muc it's betweenness centrality
//    muc_subgraph.remove_edge(src, dst);
//    tr1_map_t(int) d_src_vec, d_dst_vec;
//    muc_subgraph.find_sssp(src, d_src_vec);
//    muc_subgraph.find_sssp(dst, d_dst_vec);
//    
//    //this must be commented in general.. if it's there it makes bc_map contain deltas
//    //bcc_subgraph.i_fill_map<double>(bc_map, 0);
//    
//    for(graph_hash_t::nodes_map_t::iterator
//                        it =  muc_subgraph.nodes_map.begin();
//                        it != muc_subgraph.nodes_map.end();
//                        ++it) {
//        node_id_t s = it->first;
//        if(d_src_vec[s] != d_dst_vec[s]) {
//            i_iteration(s, 
//                    src, 
//                    dst, 
//                    d_src_vec[s], 
//                    d_dst_vec[s], 
//                    bc_map);
//        }
//    }
//    muc_subgraph.insert_edge(src, dst);
//}

void muc_t::compute_bc_inc(tr1_map_t(double)& bc_map,
                        node_id_t src,
                        node_id_t dst,
                        int max_iter)
{
    if(!valid)
        return;

    muc_fast_subgraph.fill_graph(muc_subgraph);
    tmp_conn_vertex_map.clear();
    tmp_subgraph_map.clear();
    for(conn_vert_map_t::iterator it = conn_vertex_map.begin();
            it != conn_vertex_map.end();
            ++it) {
        node_id_t new_v = muc_fast_subgraph.outin_label_map[it->first];
        tmp_conn_vertex_map.insert(make_pair(new_v, it->second));
    }
    for(subgraph_map_t::iterator it = subgraph_map.begin();
            it != subgraph_map.end();
            ++it) {
        node_id_t new_v = muc_fast_subgraph.outin_label_map[it->first];
        tmp_subgraph_map.insert(make_pair(new_v, it->second));
    }
    
    src = muc_fast_subgraph.outin_label_map[src];
    dst = muc_fast_subgraph.outin_label_map[dst];
    
    //bc_map will have for each vertex in the muc it's betweenness centrality
    muc_fast_subgraph.remove_edge(src, dst);
    vector<int> d_src_vec, d_dst_vec;
    muc_fast_subgraph.find_sssp(src, d_src_vec);
    muc_fast_subgraph.find_sssp(dst, d_dst_vec);
    
    //this must be commented in general.. if it's there it makes bc_map contain deltas
    //bcc_subgraph.i_fill_map<double>(bc_map, 0);
    if(max_iter == -1)
        max_iter = muc_fast_subgraph.nodes_vec.size();
    else
        max_iter = min(max_iter, (int)muc_fast_subgraph.nodes_vec.size());
    for(int i = 0; i < max_iter; ++i) {
        node_id_t s = i;
        if(d_src_vec[s] != d_dst_vec[s]) {
            i_iteration(s, 
                    src, 
                    dst, 
                    d_src_vec[s], 
                    d_dst_vec[s], 
                    bc_map);
        }
    }
    muc_fast_subgraph.insert_edge(src, dst);
}



void muc_t::i_iteration(node_id_t s,
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
    } else {
        //dbg_iteration(s, src, dst, d_src, d_dst, bc_vec);
        i_iteration_2(s, src, dst, d_src, d_dst, bc_vec);
    }
}


void muc_t::i_iteration_1(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_map)
{
        vector<vector<node_id_t> >&  P = muc_fast_subgraph.P; 
        vector<int>&                 path_cnt_map = muc_fast_subgraph.path_cnt_vec;
        vector<int>&                 new_path_cnt_map = muc_fast_subgraph.new_path_cnt_vec;
        vector<int>&                 dist_map = muc_fast_subgraph.dist_vec;
        vector<double>&              pair_dep_map = muc_fast_subgraph.pair_dep_vec;
        vector<double>&              new_pair_dep_map = muc_fast_subgraph.new_pair_dep_vec;
        vector<double>&              sigma_t_map = muc_fast_subgraph.sigma_t_map;
        vector<double>&              new_sigma_t_map = muc_fast_subgraph.new_sigma_t_map; 
        queue<node_id_t>&            Q = muc_fast_subgraph.Q;
        stack<node_id_t>             S;
        
        muc_fast_subgraph.init_maps();           

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
            for(int i = 0; i < muc_fast_subgraph.nodes_vec[v_i].size(); ++i) {
                node_id_t v_n = muc_fast_subgraph.nodes_vec[v_i][i];
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
            if(   tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()
               && tmp_conn_vertex_map.find(v_n) != tmp_conn_vertex_map.end()
               && s != v_n) {
                int VG_s, VG_n;
                VG_s = tmp_subgraph_map[s].size();
                VG_n = tmp_subgraph_map[v_n].size();
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
                if(tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()) {
                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
                    new_sigma_t_map[v_p] = new_sigma_t_map[v_p] + new_sigma_t_map[v_n]*new_sp_sn;
                    
                    node_id_t new_v_p = muc_fast_subgraph.inout_label_map[v_p];
                    bc_map[new_v_p] = bc_map[new_v_p] - sigma_t_map[v_n]*sp_sn/2.0;
                    bc_map[new_v_p] = bc_map[new_v_p] + new_sigma_t_map[v_n]*new_sp_sn/2.0;
                }
            }
            //IMP: this is the only change that happens to P, @src should be added as parent for dst
            if(v_n == dst) {
                node_id_t v_p = src;
                double new_sp_sn = ((double)new_path_cnt_map[v_p]/new_path_cnt_map[v_n]);
                new_pair_dep_map[v_p] = new_pair_dep_map[v_p] + new_sp_sn*(1+new_pair_dep_map[v_n]);
                if(tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()) {
                    new_sigma_t_map[v_p] = new_sigma_t_map[v_p] + new_sigma_t_map[v_n]*new_sp_sn;
                    
                    node_id_t new_v_p = muc_fast_subgraph.inout_label_map[v_p];
                    bc_map[new_v_p] = bc_map[new_v_p] + new_sigma_t_map[v_n]*new_sp_sn/2.0;
                }
            }
            if(s != v_n) {
                node_id_t new_v_n = muc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] - pair_dep_map[v_n]/2.0;
                bc_map[new_v_n] = bc_map[new_v_n] + new_pair_dep_map[v_n]/2.0;
            }
            if(tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()) {
                int VG_s = tmp_subgraph_map[s].size();
                
                node_id_t new_v_n = muc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] - pair_dep_map[v_n]*VG_s;
                bc_map[new_v_n] = bc_map[new_v_n] + new_pair_dep_map[v_n]*VG_s;
            }
        }
}

void muc_t::i_iteration_2(node_id_t s,
                        node_id_t src,
                        node_id_t dst,
                        int d_src,
                        int d_dst,
                        tr1_map_t(double)& bc_map)
{
        vector<vector<node_id_t> >&  P = muc_fast_subgraph.P; 
        vector<int>&                 path_cnt_map = muc_fast_subgraph.path_cnt_vec;
        vector<int>&                 path_cnt_inc_map = muc_fast_subgraph.path_cnt_inc_vec;
        vector<int>&                 dist_map = muc_fast_subgraph.dist_vec;
        vector<double>&              pair_dep_map = muc_fast_subgraph.pair_dep_vec;
        vector<double>&              sigma_t_map = muc_fast_subgraph.sigma_t_map;
        vector<bool>&                visited_vec = muc_fast_subgraph.visited_vec;
        queue<node_id_t>&            Q = muc_fast_subgraph.Q;
        vector<node_id_t>            S;
        
        muc_fast_subgraph.init_maps();           

        path_cnt_map[s] = 1;
        dist_map[s] = 0;
        Q.push(s);
        while(!Q.empty()) {
            node_id_t v_i = Q.front(); Q.pop();
            S.push_back(v_i);
            for(int i = 0; i < muc_fast_subgraph.nodes_vec[v_i].size(); ++i) {
                node_id_t v_n = muc_fast_subgraph.nodes_vec[v_i][i];
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
            if(   tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()
               && tmp_conn_vertex_map.find(v_n) != tmp_conn_vertex_map.end()
               && s != v_n) {
                int VG_s, VG_n;
                VG_s = tmp_subgraph_map[s].size();
                VG_n = tmp_subgraph_map[v_n].size();
                int c_t = VG_s*VG_n;
                sigma_t_map[v_n] =  sigma_t_map[v_n] + c_t;
                /*this guy must not change!*/
                //bc_map[v_n] = bc_map[v_n] - c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)path_cnt_map[v_p]/path_cnt_map[v_n]);
                pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn*(1+pair_dep_map[v_n]);
                if(tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()) {
                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
                    
                    node_id_t new_v_p = muc_fast_subgraph.inout_label_map[v_p];
                    bc_map[new_v_p] = bc_map[new_v_p] - sigma_t_map[v_n]*sp_sn/2.0;
                }
            }
            if(s != v_n) {
                node_id_t new_v_n = muc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] - pair_dep_map[v_n]/2.0;
            }
            if(tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()) {
                int VG_s = tmp_subgraph_map[s].size();
                
                node_id_t new_v_n = muc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] - pair_dep_map[v_n]*VG_s;
            }
        }
                
        //Compute new path counts and paths
        Q.push(dst);
        dist_map[dst] = dist_map[src] + 1;
        P[dst].clear();
        P[dst].push_back(src);
        path_cnt_inc_map[dst] = path_cnt_map[src];
        path_cnt_map[dst] = path_cnt_inc_map[dst];
        visited_vec[dst] = true;
        
        while (!Q.empty()) {
            node_id_t v = Q.front();
            Q.pop();
            vector<node_id_t> nbr_vec = muc_fast_subgraph.get_nbrs(v);
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
        
        muc_fast_subgraph.i_fill_map<double>(pair_dep_map, 0);
        muc_fast_subgraph.i_fill_map<double>(sigma_t_map, 0);
        //RBFS to add the new pair dependencies
        for(int i = S.size()-1; i >= 0; --i) {
            node_id_t v_n = S[i];
            if(   tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()
               && tmp_conn_vertex_map.find(v_n) != tmp_conn_vertex_map.end()
               && s != v_n) {
                int VG_s, VG_n;
                VG_s = tmp_subgraph_map[s].size();
                VG_n = tmp_subgraph_map[v_n].size();
                int c_t = VG_s*VG_n;
                sigma_t_map[v_n] =  sigma_t_map[v_n] + c_t;
                /*this guy must not change!*/
                //bc_map[v_n] = bc_map[v_n] + c_t;
            }
            for(int i = 0; i < P[v_n].size(); ++i) {
                node_id_t v_p = P[v_n][i];
                double sp_sn = ((double)path_cnt_map[v_p]/path_cnt_map[v_n]);
                pair_dep_map[v_p] = pair_dep_map[v_p] + sp_sn*(1+pair_dep_map[v_n]);
                if(tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()) {
                    sigma_t_map[v_p] = sigma_t_map[v_p] + sigma_t_map[v_n]*sp_sn;
                    
                    node_id_t new_v_p = muc_fast_subgraph.inout_label_map[v_p];
                    bc_map[new_v_p] = bc_map[new_v_p] + sigma_t_map[v_n]*sp_sn/2.0;
                }
            }
            if(s != v_n) {
                node_id_t new_v_n = muc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] + pair_dep_map[v_n]/2.0;
            }
            if(tmp_conn_vertex_map.find(s) != tmp_conn_vertex_map.end()) {
                int VG_s = tmp_subgraph_map[s].size();
                
                node_id_t new_v_n = muc_fast_subgraph.inout_label_map[v_n];
                bc_map[new_v_n] = bc_map[new_v_n] + pair_dep_map[v_n]*VG_s;
            }
        }
}