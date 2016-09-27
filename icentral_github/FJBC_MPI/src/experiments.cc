
#include <algorithm>
#include <cmath>

#include "experiments.h"
#include "stdlib.h"
#include "utility.h"
#include "bc.h"



double fuad_ideal_speedup(double iter_frac,
        int BCC_n,
        int BCC_m,
        int G_n,
        int G_m,
        int inc_bfs_ratio)
{
   double res;
   int BCC_s = BCC_n + BCC_m;
   int G_s = G_n + G_m;
   res =  1.0/(((double)BCC_n/G_n) * iter_frac * ((double)BCC_s/G_s));
   res = res * 100.0/inc_bfs_ratio; 
   return res;
}

/*
 * generates @num_edges random edges that are not in @graph
 */
void gen_rand_edges(int num_edges,
        graph_t& graph,
        vector<edge_t>& out_vec)
{
    out_vec.clear();
    for(int i = 0; i < num_edges; ++i) {
        edge_t rand_edge;
        node_id_t src, dst;
        do {
            src = rand()%graph.size();
            dst = rand()%graph.size();
            rand_edge.first    = src;
            rand_edge.second   = dst;
        } while(graph.has_edge(rand_edge) ||
                find(out_vec.begin(), out_vec.end(), rand_edge) != out_vec.end());
        out_vec.push_back(rand_edge);
    }
}

/*
 * generates @num_edges random edges that are in @graph, but are not bridges
 */
void gen_rand_edges_deletions(int num_edges,
        graph_t& graph,
        vector<edge_t>& out_vec)
{
    vector<edge_t> bridges_vec;
    graph.find_bridge_edges(bridges_vec);
    //NOTE:insert each edge twice -- e(v1, v2) e(v2, v1)
    int num_bridges = bridges_vec.size();
    for(int i = 0; i < num_bridges; ++i) {
        bridges_vec.push_back(make_pair(bridges_vec[i].second, bridges_vec[i].first));
    }
    out_vec.clear();
    int num_edges_graph = graph.edge_set.size();
    set<edge_t>::const_iterator it(graph.edge_set.begin());
    for(int i = 0; i < num_edges; ++i) {
        edge_t rand_edge;
        do {
            //select an edge at random
            it = graph.edge_set.begin();
            advance(it, rand()%num_edges_graph);
            rand_edge = *it;
        } while(find(bridges_vec.begin(), bridges_vec.end(), rand_edge) != bridges_vec.end());
        out_vec.push_back(rand_edge);
    }
}

void insertion_test_qube(graph_t& graph, int num_edges_to_insert)
{
    //1. prepare structures to store the results for each edge 
    //2. generate random edges
    //3. run Brandes and find time
    //4. run update after each edge insertion and record results
    //5. print the results progressively
    //int num_edges_to_insert = 2;
    
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    double brandes_time, muc_find_time;
    timer tm;
    
    printf("Brandes runtime: \t");
    tm.start();
    vector<double> bc_vec = brandes_betweenness_centrality(graph);
    tm.stop();
    brandes_time = tm.interval();
    printf("%f\n", brandes_time);
    
    printf("# MUCs: \t\t");
    tm.start();
    graph.find_mucs();
    tm.stop();
    muc_find_time = tm.interval();
    printf("%d\n", graph.get_num_mucs());
    printf("Time to find MUCs: \t%f\n", muc_find_time);
    
    vector<edge_t> rand_edge_vec;
    printf("Generating edges..n\n");
    for(int i = 0; i < num_edges_to_insert; ++i) {
        edge_t rand_edge;
        node_id_t src, dst;
        do {
            src = rand()%graph.size();
            dst = rand()%graph.size();
            rand_edge.first    = src;
            rand_edge.second   = dst;
        } while(graph.has_edge(rand_edge));
        rand_edge_vec.push_back(rand_edge);
    }

    vector<edge_stat_t> edge_stat_vec;
    
    //==========================================================
//    node_id_t nd = rand_edge_vec[0].first;
//    muc_id_t muc_id = graph.node_to_muc_vec[nd];
//    graph.muc_vec[muc_id].muc_subgraph.print_graph(false);
//    return;
    //==========================================================
    
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        edge_t rand_edge = rand_edge_vec[i];
        edge_stat_t edge_stat;
        printf("Adding edge (%d, %d)\n", rand_edge.first, rand_edge.second);
        tm.start();
        graph.insert_edge_update_muc(rand_edge.first, rand_edge.second);
        tm.stop();
        edge_stat.edge_ins_time = tm.interval();
        printf("\tTime to update MUCs: \t%f\n", edge_stat.edge_ins_time);

        muc_id_t muc_id = graph.node_to_muc_vec[rand_edge.first];
        edge_stat.muc_num_edges = graph.muc_vec[muc_id].muc_subgraph.edge_set.size();
        edge_stat.muc_num_nodes = graph.muc_vec[muc_id].muc_subgraph.nodes_map.size();
        tr1::unordered_map<node_id_t, double> bc_map;
        tm.start();
        graph.muc_vec[muc_id].compute_bc(bc_map);
        tm.stop();
        edge_stat.muc_bc_update_time = tm.interval();
        printf("\tTime to update BC: \t%f\n", edge_stat.muc_bc_update_time);
        printf("\tSpeedup: \t\t%f\n", brandes_time/edge_stat.muc_bc_update_time);
        printf("\tMUC # nodes: \t\t%d\n", edge_stat.muc_num_nodes);
        printf("\tMUC # edges: \t\t%d\n", edge_stat.muc_num_edges);
        edge_stat_vec.push_back(edge_stat);
    }
    
    
    double total_muc_update_time = 0;
    for(int i = 0; i < edge_stat_vec.size(); ++i) {
        total_muc_update_time += edge_stat_vec[i].muc_bc_update_time;
    }
    double avg_speedup = (num_edges_to_insert*brandes_time)/total_muc_update_time;
    
    printf("[[Average speedup: %f]]\n", avg_speedup);
}


//TODO change name.. not hash anymore
void insertion_test_qube_hash(graph_t& graph,
        int num_edges_to_insert,
        vector<edge_t> rand_edge_vec,
        bool use_incremental_algo,
        int max_iter)
{
    //1. prepare structures to store the results for each edge 
    //2. generate random edges
    //3. run Brandes and find time
    //4. run update after each edge insertion and record results
    //5. print the results progressively
    //int num_edges_to_insert = 2;
    
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    graph_hash_t g;
    g.fill_graph(graph);
    
    double brandes_time, muc_find_time;
    timer tm;
    
    printf("Brandes runtime: \t");
    tr1_map_t(double) bc_vec;
    tm.start();
    bc_vec = brandes_bc_hash_out(graph, max_iter);
    tm.stop();
    if(max_iter == -1)
        brandes_time = tm.interval();
    else
        brandes_time = (tm.interval()/max_iter)*graph.size();
    printf("%f\n", brandes_time);
    printf("Brandes avg iter time: \t%f\n", brandes_time/graph.size());
    
    printf("# MUCs: \t\t");
    tm.start();
    graph.find_mucs();
    tm.stop();
    muc_find_time = tm.interval();
    printf("%d\n", graph.get_num_mucs());
    printf("Time to find MUCs: \t%f\n", muc_find_time);
    
    //vector<edge_t> rand_edge_vec;
    if(rand_edge_vec.size() == 0) {
        printf("Generating edges...\n");
        gen_rand_edges(num_edges_to_insert, graph, rand_edge_vec);
    }
    vector<edge_stat_t> edge_stat_vec;
    
    double tot_speedup = 0;
    double tot_ideal_speedup = 0;
    
    vector<double> speedup_vec;
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        edge_t rand_edge = rand_edge_vec[i];
        edge_stat_t edge_stat;
        printf("Adding edge (%d, %d)\n", rand_edge.first, rand_edge.second);
        tm.start();
        graph.insert_edge_update_muc(rand_edge.first, rand_edge.second);
        tm.stop();
        edge_stat.edge_ins_time = tm.interval();

        muc_id_t muc_id = graph.node_to_muc_vec[rand_edge.first];
        edge_stat.muc_num_edges = graph.muc_vec[muc_id].muc_subgraph.edge_set.size();
        edge_stat.muc_num_nodes = graph.muc_vec[muc_id].muc_subgraph.nodes_map.size();
        tr1::unordered_map<node_id_t, double> bc_map;
        tm.start();
        if(!use_incremental_algo) {
            graph.muc_vec[muc_id].compute_bc(bc_map, max_iter);
        } else {
            graph.muc_vec[muc_id].compute_bc_inc(bc_map, rand_edge.first, rand_edge.second);
        }
        tm.stop();
        if(max_iter == -1)
            edge_stat.muc_bc_update_time = tm.interval();
        else
            edge_stat.muc_bc_update_time = (tm.interval()/max_iter)*edge_stat.muc_num_nodes;
        //TMP CODE XXX FIX THIS SHIT
        edge_stat.muc_bc_update_time += edge_stat.edge_ins_time;
        ////
        
        double speedup = brandes_time/edge_stat.muc_bc_update_time;
        tot_speedup += speedup;
        speedup_vec.push_back(speedup);
        //estimating ideal speedup using the node and edge proportions
        double ideal_speedup;
        double V  = graph.size(); double E = graph.edge_set.size();
        double Vm = edge_stat.muc_num_nodes; double Em = edge_stat.muc_num_edges;
        ideal_speedup = ((V+E)/(Vm+Em))*(V/Vm);
        tot_ideal_speedup += ideal_speedup;
        ////
        
        printf("\tMUC # nodes: \t\t%d (%f G_n)\n",
                edge_stat.muc_num_nodes,
                edge_stat.muc_num_nodes/(double)graph.size());
        printf("\tMUC # edges: \t\t%d (%f G_M)\n",
                edge_stat.muc_num_edges,
                edge_stat.muc_num_edges/(double)graph.edge_set.size());
        printf("\n");
        printf("\tTime to update BC: \t%f\n", edge_stat.muc_bc_update_time);
        printf("\tSpeedup: \t\t%f (Ideal: %f)\n", speedup, ideal_speedup);
        printf("\n");
        printf("\tQUBE avg iter time: \t%f\n",
                edge_stat.muc_bc_update_time/edge_stat.muc_num_nodes);
        printf("\tTime to update MUCs: \t%f\n", edge_stat.edge_ins_time);

        edge_stat_vec.push_back(edge_stat);
    }
    
    double avg_speedup = tot_speedup/num_edges_to_insert;
    double avg_ideal_speedup = tot_ideal_speedup/num_edges_to_insert;
    printf("[[Average speedup: %f (Ideal: %f)]]\n",
        avg_speedup,
        avg_ideal_speedup);
    
    //calculating 95% confidence interval
    double s_mean, s_median, s_stddev, s_min, s_max;
    simple_stats(speedup_vec, s_mean, s_stddev, s_median, s_min, s_max);
    double diff = 1.96*s_stddev/sqrt(num_edges_to_insert);
    double v1 = s_mean - diff;
    double v2 = s_mean + diff;
    printf("[[Speedup 95 CI: %f - %f]]\n", v1, v2);
}


void incremental_brandes_test(graph_t& graph, int num_edges_to_insert, vector<edge_t> rand_edge_vec)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    double brandes_time;
    timer tm;
    
    printf("Brandes avg time per iteration: \t");
    tm.start();
    vector<double> bc_vec = brandes_bc(graph);
    tm.stop();
    brandes_time = tm.interval();
    printf("%f\n", brandes_time/graph.size());
    printf("Brandes time: \t\t\t\t%f\n", brandes_time);
    
    if(rand_edge_vec.size() == 0) {
        printf("Generating edges...\n");
        gen_rand_edges(num_edges_to_insert, graph, rand_edge_vec);
    }

    vector<edge_stat_t> edge_stat_vec;
 
    vector<vector<int> >    cnt_vec2;
    vector<vector<double> > time_vec2;
    
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        //printf("Inserting edge [%d]...\n", i);
        edge_t rand_edge = rand_edge_vec[i];
        
        vector<int> cnt_vec;
        vector<double> time_vec;
        cnt_vec.resize(3);
        time_vec.resize(3);
        
        
        incremental_brandes_experimental(graph,
                rand_edge.first, 
                rand_edge.second, 
                bc_vec,
                time_vec,
                cnt_vec);
        cnt_vec2.push_back(cnt_vec);
        time_vec2.push_back(time_vec);
    }

    int tot_cnt_arr[] = {0, 0, 0};
    double tot_time_arr[] = {0, 0, 0};
    double tot_avg_arr[] = {0, 0, 0};
    
    for(int i = 0; i < cnt_vec2.size(); ++i) {
        tot_cnt_arr[0] += cnt_vec2[i][0];
        tot_cnt_arr[1] += cnt_vec2[i][1];
        tot_cnt_arr[2] += cnt_vec2[i][2];
        
        tot_time_arr[0] += time_vec2[i][0];
        tot_time_arr[1] += time_vec2[i][1];
        tot_time_arr[2] += time_vec2[i][2];
        
        tot_avg_arr[0] += (time_vec2[i][0]/cnt_vec2[i][0]);
        tot_avg_arr[1] += (time_vec2[i][1]/cnt_vec2[i][1]);
        tot_avg_arr[2] += (time_vec2[i][2]/cnt_vec2[i][2]);
    }
    double inc_tot_time = tot_time_arr[0] + tot_time_arr[1] + tot_time_arr[2];
    double inc_avg_time = inc_tot_time/rand_edge_vec.size();
    printf("Incremental Brandes time:\t\t%f\n", inc_avg_time);
    printf("\n");
    double avg_d1_tm = tot_time_arr[1]/tot_cnt_arr[1];
    double avg_d2_tm = tot_time_arr[2]/tot_cnt_arr[2];
    double avg_br_tm = brandes_time/graph.size();
    printf("\n");
    printf("Tot d0 BFS:\t\t%d\n", tot_cnt_arr[0]);
    printf("Tot d1 BFS:\t\t%d\n", tot_cnt_arr[1]);
    printf("Tot d2 BFS:\t\t%d\n", tot_cnt_arr[2]);
    printf("\n");
    printf("Avg d1 BFS time:\t\t%f (%f)\n", avg_d1_tm, avg_d1_tm/avg_br_tm);
    printf("Avg d2 BFS time:\t\t%f (%f)\n", avg_d2_tm, avg_d2_tm/avg_br_tm);
    printf("Avg Brandes time:\t\t%f\n", avg_br_tm);
    double avg_inc_iter = 
    (tot_time_arr[1] + tot_time_arr[2])/(tot_cnt_arr[1] + tot_cnt_arr[2]);
    
    double mem_dag_time = (avg_d1_tm-avg_br_tm)*tot_cnt_arr[1] + (avg_d2_tm-avg_br_tm)*tot_cnt_arr[2];
    double mem_dag_speedup = brandes_time/mem_dag_time;
    printf("Ideal mem DAG speedup:\t %f\n", mem_dag_speedup);
    printf("Avg inc BFS time:\t\t%f (%f)\n", avg_inc_iter, avg_inc_iter/avg_br_tm);
    printf("[[Avg speedup:\t\t\t%f]]\n", brandes_time/inc_avg_time);
    
}

void insertion_test_fuad_hash(graph_t& graph, int num_edges_to_insert)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    graph_hash_t g;
    g.fill_graph(graph);
    
    double brandes_time;
    timer tm;
    
    int num_iter = g.size();
    printf("Brandes runtime: \t");
    tm.start();
    tr1_map_t(double) bc_vec = brandes_bc_hash(g);
    tm.stop();
    brandes_time = tm.interval();
    printf("%f\n", brandes_time);
    printf("Brandes avg iter time: \t%f\n", brandes_time/num_iter);
    
    vector<edge_t> rand_edge_vec;
    printf("Generating edges...\n");
    for(int i = 0; i < num_edges_to_insert; ++i) {
        edge_t rand_edge;
        node_id_t src, dst;
        do {
            src = rand()%graph.size();
            dst = rand()%graph.size();
            rand_edge.first    = src;
            rand_edge.second   = dst;
        } while(graph.has_edge(rand_edge) ||
                find(rand_edge_vec.begin(), rand_edge_vec.end(), rand_edge) != rand_edge_vec.end());
        rand_edge_vec.push_back(rand_edge);
    }

    double tot_fuad_time = 0;
    double tot_speedup = 0;
    
    vector<bcc_stat_t> bcc_stat_vec;
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        bcc_stat_t bcc_stat;
        edge_t rand_edge = rand_edge_vec[i];
        printf("Inserting edge (%d, %d) and updating bc...\n", rand_edge.first, rand_edge.second);
        tm.start();
        graph.insert_edge_update_bc_experimental(rand_edge.first, rand_edge.second, bcc_stat);
        tm.stop();
        
        double avg_iter_time_naive;
        graph.approx_bcc_iter_tm(rand_edge.first, rand_edge.second, avg_iter_time_naive);
        tot_fuad_time += tm.interval();
        
        tot_speedup += (brandes_time/tm.interval());
        
        printf("\tBCC num nodes: \t\t%d\n", bcc_stat.bcc_num_nodes);
        printf("\tBCC num edges: \t\t%d\n", bcc_stat.bcc_num_edges);
        printf("\n");
        printf("\tTime to update BC: \t%f\n", tm.interval());
        printf("\tSpeedup: \t\t%f\n", brandes_time/tm.interval());
        printf("\n");
        printf("\tBCC num d0 BFS: \t%d\n", bcc_stat.num_d0_iter);
        printf("\tBCC d0 time: \t\t%f\n", bcc_stat.tot_d0_tm/bcc_stat.num_d0_iter);
        printf("\tBCC num d1 BFS: \t%d\n", bcc_stat.num_d1_iter);
        printf("\tBCC d1 time: \t\t%f\n", bcc_stat.tot_d1_tm/bcc_stat.num_d1_iter);
        printf("\tBCC num d2 BFS: \t%d\n", bcc_stat.num_d2_iter);
        printf("\tBCC d2 time: \t\t%f\n", bcc_stat.tot_d2_tm/bcc_stat.num_d2_iter);
        printf("\tavg iter time naive: \t%f\n", avg_iter_time_naive);
        //printf("\tBCC bcc find time: \t%f\n", bcc_stat.bcc_find_time);
        //printf("\tBCC sssp time: \t\t%f\n", bcc_stat.sssp_tm);
        
        
        
        bcc_stat_vec.push_back(bcc_stat);
    }
    
    //double avg_speedup = (num_edges_to_insert*brandes_time)/tot_fuad_time;
    double avg_speedup = tot_speedup/num_edges_to_insert;
    
    printf("[[Average speedup: %f]]\n", avg_speedup); 
}

void insertion_test_fuad_max_iter(graph_t& graph,
        int num_edges_to_insert,
        int num_iter,
        double time_to_run_all_iter,
        vector<edge_t> rand_edge_vec)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    
    double brandes_time;
    timer tm;
    
    printf("Brandes runtime: \t");
    tm.start();
    tr1_map_t(double) bc_vec = brandes_bc_hash_out(graph, num_iter);
    tm.stop();
    brandes_time = (tm.interval()*graph.size()/num_iter);
    
    bool do_all_iter = false;
    if(brandes_time < time_to_run_all_iter)
        do_all_iter = true;
    
    if(do_all_iter) {
        num_iter = -1;
        tm.start();
        tr1_map_t(double) bc_vec = brandes_bc_hash_out(graph, num_iter);
        tm.stop();
        brandes_time = tm.interval();
        printf("%f\n", brandes_time);
        printf("Brandes avg iter time: \t%f\n", tm.interval()/graph.size());
        printf("WILL DO ALL ITERATIONS...\n");
    } else {
        printf("%f\n", brandes_time);
        printf("Brandes avg iter time: \t%f\n", tm.interval()/num_iter);
    }

    
    //vector<edge_t> rand_edge_vec;
    if(rand_edge_vec.size() == 0) {
        printf("Generating edges...\n");
        gen_rand_edges(num_edges_to_insert, graph, rand_edge_vec);
    }

    double tot_fuad_time = 0;
    double tot_speedup = 0;
    
    double tot_ideal_speedup_perfect = 0;
    double tot_ideal_speedup = 0;
    
    vector<bcc_stat_t>  bcc_stat_vec;
    vector<double>      speedup_vec;
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        bcc_stat_t bcc_stat;
        edge_t rand_edge = rand_edge_vec[i];
        double est_time = 0;
        
        printf("Inserting edge (%d, %d) and updating bc...\n", rand_edge.first, rand_edge.second);
        tm.start();
        graph.insert_edge_update_bc_experimental(rand_edge.first, rand_edge.second, bcc_stat, num_iter, num_iter);
        tm.stop();
        if(do_all_iter) {
            est_time = bcc_stat.bc_update_time + bcc_stat.bcc_find_time;
            
            // TMP XXX FIX THIS
            //est_time = bcc_stat.bc_update_time;
            ////
        }
        
        double avg_iter_time_naive;
        graph.approx_bcc_iter_tm(rand_edge.first, rand_edge.second, avg_iter_time_naive, num_iter);
        tot_fuad_time += tm.interval();
        
        //fix counts if will produce errors:
        if(bcc_stat.num_d0_iter == 0) {
            bcc_stat.num_d0_iter = 1;
            bcc_stat.tot_d0_tm = 0;
        }
        if(bcc_stat.num_d1_iter == 0) {
            bcc_stat.num_d1_iter = 1;
            bcc_stat.tot_d1_tm = 0;
        }
        if(bcc_stat.num_d2_iter == 0) {
            bcc_stat.num_d2_iter = 1;
            bcc_stat.tot_d2_tm = 0;
        }
        ////
        if(!do_all_iter) {
            est_time += bcc_stat.tot_d1_tm*bcc_stat.tot_d1_iter/bcc_stat.num_d1_iter;
            est_time += bcc_stat.tot_d2_tm*bcc_stat.tot_d2_iter/bcc_stat.num_d2_iter;
        }
        

        tot_speedup += (brandes_time/est_time);
        speedup_vec.push_back(brandes_time/est_time);
        
        //ideal speedup estimation
        double ideal_speedup_perfect, ideal_speedup;
        double G_n   = graph.size(); 
        double G_m   = graph.edge_set.size();
        double BCC_n = bcc_stat.bcc_num_nodes;
        double BCC_m = bcc_stat.bcc_num_edges;
        double inc_iter_ratio = 150;//PARAMETER!! IMP make a #define
        double iter_frac = 1 - ((double)bcc_stat.tot_d0_iter/BCC_n);
        ideal_speedup_perfect =
                fuad_ideal_speedup(iter_frac, 
                    BCC_n,
                    BCC_m,
                    G_n,
                    G_m,
                    100);
        ideal_speedup =
                fuad_ideal_speedup(iter_frac, 
                    BCC_n,
                    BCC_m,
                    G_n,
                    G_m,
                    inc_iter_ratio);
        
        tot_ideal_speedup_perfect += ideal_speedup_perfect;
        tot_ideal_speedup += ideal_speedup;
        ////
        
        //speedup if BFS DAGS were to be stored in memory
        double in_mem_speedup;
        double in_mem_time = 0;
        in_mem_time += (bcc_stat.tot_d1_tm/bcc_stat.num_d1_iter-avg_iter_time_naive)*bcc_stat.tot_d1_iter;
        in_mem_time += (bcc_stat.tot_d2_tm/bcc_stat.num_d2_iter-avg_iter_time_naive)*bcc_stat.tot_d2_iter;
        in_mem_speedup = brandes_time/in_mem_time;
        printf("\tStored DAGS ideal speedup: \t%f\n\n", in_mem_speedup);
        ////
        
        printf("\tBCC num nodes: \t\t\t%d (%f G_n)\n",
                bcc_stat.bcc_num_nodes,
                (double)bcc_stat.bcc_num_nodes/graph.size());
        printf("\tBCC num edges: \t\t\t%d (%f G_m)\n",
                bcc_stat.bcc_num_edges,
                (double)bcc_stat.bcc_num_edges/graph.edge_set.size());
        printf("\n");
        printf("\tTime to update BC: \t\t%f\n",
                est_time);
        printf("\tSpeedup: \t\t\t%f (Ideal: %f - %f)\n",
                brandes_time/est_time,
                ideal_speedup,
                ideal_speedup_perfect);
        printf("\n");
        printf("\tBCC d1 time: \t\t\t%f (%f to naive iter)\n",
                bcc_stat.tot_d1_tm/bcc_stat.num_d1_iter,
                (bcc_stat.tot_d1_tm/bcc_stat.num_d1_iter)/avg_iter_time_naive);
        printf("\tBCC d2 time: \t\t\t%f (%f to naive iter)\n",
                bcc_stat.tot_d2_tm/bcc_stat.num_d2_iter,
                (bcc_stat.tot_d2_tm/bcc_stat.num_d2_iter)/avg_iter_time_naive);
        printf("\tavg iter time naive: \t\t%f (%f to naive iter)\n",
                avg_iter_time_naive,
                avg_iter_time_naive/avg_iter_time_naive);
        printf("\n");
        printf("\tBCC num/tot d0 BFS: \t\t%*d/%*d (%f to BCC_n)\n",
                8,
                bcc_stat.num_d0_iter,
                8,
                bcc_stat.tot_d0_iter,
                bcc_stat.tot_d0_iter/(double)bcc_stat.bcc_num_nodes);
        printf("\tBCC num/tot d1 BFS: \t\t%*d/%*d (%f to BCC_n)\n",
                8,
                bcc_stat.num_d1_iter,
                8,
                bcc_stat.tot_d1_iter,
                bcc_stat.tot_d1_iter/(double)bcc_stat.bcc_num_nodes);
        printf("\tBCC num/tot d2 BFS: \t\t%*d/%*d (%f to BCC_n)\n",
                8,
                bcc_stat.num_d2_iter,
                8,
                bcc_stat.tot_d2_iter,
                bcc_stat.tot_d2_iter/(double)bcc_stat.bcc_num_nodes);
        printf("\n");
        printf("\tBCC bcc find time: \t\t%f\n", bcc_stat.bcc_find_time);
        printf("\tBCC sssp time: \t\t\t%f\n", bcc_stat.sssp_tm);
        
        bcc_stat_vec.push_back(bcc_stat);
        
        //dont remove the edge, to have fair comparison to QUBE
        //in my QUBE implementation the MUC update after edge removal is not
        //implemented
        //graph.remove_edge(rand_edge.first, rand_edge.second);
    }
    
    //double avg_speedup = (num_edges_to_insert*brandes_time)/tot_fuad_time;
    double avg_speedup = tot_speedup/num_edges_to_insert;
    double avg_ideal_speedup_perfect = tot_ideal_speedup_perfect/num_edges_to_insert;
    double avg_ideal_speedup = tot_ideal_speedup/num_edges_to_insert;
    
    printf("[[Average speedup: %f (Ideal: %f - %f)]]\n",
        avg_speedup,
        avg_ideal_speedup,
        avg_ideal_speedup_perfect);
    
    //calculating 95% confidence interval
    double s_mean, s_median, s_stddev, s_min, s_max;
    simple_stats(speedup_vec, s_mean, s_stddev, s_median, s_min, s_max);
    double diff = 1.96*s_stddev/sqrt(num_edges_to_insert);
    double v1 = s_mean - diff;
    double v2 = s_mean + diff;
    printf("[[Speedup 95 CI: %f - %f]]\n", v1, v2);
}


void compare_brandes_hash_vs_vector(graph_t& graph, int num_iter)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    graph_hash_t g;
    g.fill_graph(graph);
    
    double brandes_time;
    timer tm;

    printf("\n");
    printf("Brandes runtime: \t\t");
    tm.start();
    vector<double> bc_vec = brandes_bc(graph, num_iter);
    tm.stop();
    brandes_time = (tm.interval()*g.size()/num_iter);
    printf("%f\n", brandes_time);
    printf("Brandes avg iter time: \t\t%f\n", tm.interval()/num_iter);
    printf("\n");
    
    printf("Brandes hash runtime: \t\t");
    tm.start();
    tr1_map_t(double) bc_map = brandes_bc_hash(g, num_iter);
    tm.stop();
    brandes_time = (tm.interval()*g.size()/num_iter);
    printf("%f\n", brandes_time);
    printf("Brandes hash avg iter time: \t%f\n", tm.interval()/num_iter);
}

/*
 * estimates speedups assuming all edges are going to be picked in the largest
 * biconnected component
 */
void speedup_info_lbcc(graph_t& graph, int num_edges_to_insert)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
//    graph_hash_t g;
//    g.fill_graph(graph);
    vector<graph_hash_t> bcc_vec;
    graph.find_bicon_comp(bcc_vec);
    
    //find the bcc with the largest number of nodes
    int lbcc_i = 0;
    int lbcc_s = bcc_vec[0].size() + bcc_vec[0].edge_set.size();
    if(bcc_vec.size() > 1) {
        for(int i = 1; i < bcc_vec.size(); ++i) {
            if(bcc_vec[i].size() + bcc_vec[i].edge_set.size() > lbcc_s) {
                lbcc_s = bcc_vec[i].size() + bcc_vec[i].edge_set.size();
                lbcc_i = i;
            }
        }
    }
    
    int G_s = graph.size() + graph.edge_set.size();
    printf("\n");
    printf("LBCC # nodes:\t\t%d\n", bcc_vec[lbcc_i].size());
    printf("LBCC # edges:\t\t%d\n", bcc_vec[lbcc_i].edge_set.size());
    printf("\n");
    printf("LBCC size fraction:\t%f\n", (double)lbcc_s/G_s);
    printf("LBCC node fraction:\t%f\n", (double)bcc_vec[lbcc_i].size()/graph.size());
    printf("\n");
    
    vector<edge_t> rand_edge_vec;
    printf("Generating edges...\n");
    for(int i = 0; i < num_edges_to_insert; ++i) {
        edge_t rand_edge;
        node_id_t src, dst;
        do {
            src = rand()%graph.size();
            dst = rand()%graph.size();
            rand_edge.first    = src;
            rand_edge.second   = dst;
        } while(!bcc_vec[lbcc_i].has_edge(rand_edge.first, rand_edge.second) ||
                find(rand_edge_vec.begin(), rand_edge_vec.end(), rand_edge) != rand_edge_vec.end());
        rand_edge_vec.push_back(rand_edge);
    }
    
    printf("\n");
    vector<int> d0_vec, d1_vec, d2_vec;
    vector<double> d0_frac_vec;
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        node_id_t src, dst;
        src = rand_edge_vec[i].first;
        dst = rand_edge_vec[i].second;
        int d0, d1, d2;
        bcc_vec[lbcc_i].find_pruning_counts_exp(src, dst, d0, d1, d2);
        d0_vec.push_back(d0);
        d1_vec.push_back(d1);
        d2_vec.push_back(d2);
        d0_frac_vec.push_back((double)d0/bcc_vec[lbcc_i].size());
        
        printf("d0:[%*d][%f]   d1:[%*d]    d2:[%*d] -- Edge(%d, %d)\n",
                8,
                d0, (double)d0/bcc_vec[lbcc_i].size(),
                8,
                d1,
                8,
                d2,
                rand_edge_vec[i].first,
                rand_edge_vec[i].second);
    }
    printf("\n");
    
    double min, max, median;
    double mean, stddev;
    simple_stats(d0_frac_vec, mean, stddev, median, min, max);
    printf("\n");
    printf("Mean[%f] Stddev[%f] Median[%f] Min[%f] Max[%f]\n",
        mean, stddev, median, min, max);
    
    
    int lbcc_n                 = bcc_vec[lbcc_i].size();
    int G_n                    = graph.size();
    double iter_frac_median    = 1.0 - median;
    double iter_frac_mean      = 1.0 - mean;
    double iter_frac_min       = 1.0 - min;
    double iter_frac_max       = 1.0 - max;
    double iter_frac_stddev    = 1.0 - stddev;
    
    
    double est_median_speedup = 1.0/(((double)lbcc_n/G_n) * iter_frac_median * ((double)lbcc_s/G_s));
    double est_mean_speedup   = 1.0/(((double)lbcc_n/G_n) * iter_frac_mean * ((double)lbcc_s/G_s));
    double est_min_speedup    = 1.0/(((double)lbcc_n/G_n) * iter_frac_min * ((double)lbcc_s/G_s));
    double est_max_speedup    = 1.0/(((double)lbcc_n/G_n) * iter_frac_max * ((double)lbcc_s/G_s));
    double est_stddev_speedup = 1.0/(((double)lbcc_n/G_n) * iter_frac_stddev * ((double)lbcc_s/G_s));
    
    //estimate of the time it takes to do my long delta iteration
    //to the time it takes the normal iteration
    double long_iter_perc;
    long_iter_perc = 120.0;
    printf("\n");
    printf("Estimate ideal speedups: (delta iter to orig iter [%f])\n",
        long_iter_perc);
    printf("Mean[%f] Stddev[%f] Median[%f] Min[%f] Max[%f]\n",
        est_mean_speedup * (100/long_iter_perc),
        est_stddev_speedup * (100/long_iter_perc),
        est_median_speedup * (100/long_iter_perc),
        est_min_speedup * (100/long_iter_perc),
        est_max_speedup * (100/long_iter_perc));
    
    long_iter_perc = 100.0;
    printf("Estimate ideal speedups: (delta iter to orig iter [%f])\n",
        long_iter_perc);
    printf("Mean[%f] Stddev[%f] Median[%f] Min[%f] Max[%f]\n",
        est_mean_speedup * (100/long_iter_perc),
        est_stddev_speedup * (100/long_iter_perc),
        est_median_speedup * (100/long_iter_perc),
        est_min_speedup * (100/long_iter_perc),
        est_max_speedup * (100/long_iter_perc));
    
}


/*
 * estimates speedups WITHOUT assuming that all edges are going to be picked from
 * the same biconnected component
 */
void speedup_info(graph_t& graph, int num_edges_to_insert)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    vector<graph_hash_t> bcc_vec;
    graph.find_bicon_comp(bcc_vec);
    
    int G_s = graph.size() + graph.edge_set.size();
    int G_n                    = graph.size();
    
    vector<edge_t> rand_edge_vec;
    printf("Generating edges...\n");
    for(int i = 0; i < num_edges_to_insert; ++i) {
        edge_t rand_edge;
        node_id_t src, dst;
        do {
            src = rand()%graph.size();
            dst = rand()%graph.size();
            rand_edge.first    = src;
            rand_edge.second   = dst;
        } while(graph.has_edge(rand_edge) ||
                find(rand_edge_vec.begin(), rand_edge_vec.end(), rand_edge) != rand_edge_vec.end());
        rand_edge_vec.push_back(rand_edge);
    }
    
    
    printf("\n");
    vector<double> est_speedup_vec;
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        node_id_t src = rand_edge_vec[i].first;
        node_id_t dst = rand_edge_vec[i].second;
        graph_hash_t bcc;
        graph.insert_edge(src, dst);
        graph.find_edge_bcc_subgraph(bcc, src, dst);
        graph.remove_edge(src, dst);
        
        int bcc_s = bcc.size() + bcc.edge_set.size();
        int bcc_n = bcc.size();
        
        printf("BCC#nodes:[%*d]   BCC#edges:[%*d]   BCCSizeFrac:[%f]   BCCNodeFrac:[%f]\n",
                8,
                bcc_n,
                8,
                bcc.edge_set.size(),
                (double)bcc_s/G_s,
                (double)bcc_n/graph.size());
        
        int d0, d1, d2;
        bcc.find_pruning_counts_exp(src, dst, d0, d1, d2);
        
        double iter_frac = 1.0 - ((double)d0/bcc_n);
        double est_speedup = 1.0/(((double)bcc_n/G_n) * iter_frac * ((double)bcc_s/G_s));
        
        est_speedup_vec.push_back(est_speedup);
        
        printf("d0:[%*d][%f]   d1:[%*d]    d2:[%*d] -- Edge(%d, %d)\n",
                8,
                d0,
                (double)d0/bcc_n,
                8,
                d1,
                8,
                d2,
                src,
                dst);
        
        printf("\n");
    }
    
    double min, max, median;
    double mean, stddev;
    simple_stats(est_speedup_vec, mean, stddev, median, min, max);
    //estimate of the time it takes to do my long delta iteration
    //to the time it takes the normal iteration
    double long_iter_perc;
    long_iter_perc = 150.0;
    printf("\n");
    printf("Estimate ideal speedups: (delta iter to orig iter [%f])\n",
        long_iter_perc);
    printf("Mean[%f] Stddev[%f] Median[%f] Min[%f] Max[%f]\n",
        mean * (100/long_iter_perc),
        stddev * (100/long_iter_perc),
        median * (100/long_iter_perc),
        min * (100/long_iter_perc),
        max * (100/long_iter_perc));
    
    double worst_mean = mean * (100/long_iter_perc);
    long_iter_perc = 100.0;
    printf("Estimate ideal speedups: (delta iter to orig iter [%f])\n",
        long_iter_perc);
    printf("Mean[%f] Stddev[%f] Median[%f] Min[%f] Max[%f]\n",
        mean * (100/long_iter_perc),
        stddev * (100/long_iter_perc),
        median * (100/long_iter_perc),
        min * (100/long_iter_perc),
        max * (100/long_iter_perc));
    double best_mean = mean * (100/long_iter_perc);
    printf("\n");
    printf("[[Average speedup (ideal): %f - %f]]\n",
        worst_mean,
        best_mean);
}


/*
 * counts the number of biconnected components in the graph that have 
 * at least @num_edges edges
 */
void count_bcc(graph_t& graph, int num_edges)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    vector<graph_hash_t> bcc_vec;
    graph.find_bicon_comp(bcc_vec);
    
    vector<int> bcc_sz_vec;
    int num_bcc = 0;
    if(bcc_vec.size() > 1) {
        for(int i = 0; i < bcc_vec.size(); ++i) {
            bcc_sz_vec.push_back(bcc_vec[i].edge_set.size());
            if(bcc_vec[i].edge_set.size() >= num_edges) {
                num_bcc++;
            }
        }
    }
    
    printf("\n");
    if(bcc_sz_vec.size() > 2) {
        sort(bcc_sz_vec.begin(), bcc_sz_vec.end());
        printf("Num edges in LBCC:\t[%d]\n", bcc_sz_vec.back());
        printf("Num edges in 2nd largest BCC:\t[%d]\n",
                bcc_sz_vec[bcc_sz_vec.size()-2]);
    }
    printf("\n");
    double frac_edges_in_lbcc = bcc_sz_vec.back()/(double)graph.edge_set.size();
    printf("Frac edges in LBCC:\t\t[%f]\n", frac_edges_in_lbcc);
    printf("\n");
    printf("Num BCC:\t\t\t[%d]\n", bcc_vec.size());
    printf("Num BCC with # edges >= %d:\t[%d]\n", num_edges, num_bcc);
}


void qube_ideal_speedup(graph_t& graph,
        int num_edges_to_insert,
        vector<edge_t> rand_edge_vec)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    graph_hash_t g;
    g.fill_graph(graph);
    
    double brandes_time, muc_find_time;
    timer tm;
  
    printf("# MUCs: \t\t");
    tm.start();
    graph.find_mucs();
    tm.stop();
    muc_find_time = tm.interval();
    printf("%d\n", graph.get_num_mucs());
    printf("Time to find MUCs: \t%f\n", muc_find_time);
    
    printf("Generating edges...\n");
    if(rand_edge_vec.size() == 0) {
        gen_rand_edges(num_edges_to_insert, graph, rand_edge_vec);
    }
    vector<edge_stat_t> edge_stat_vec;
    
    double tot_ideal_speedup = 0;
    
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        edge_t rand_edge = rand_edge_vec[i];
        edge_stat_t edge_stat;
        printf("Adding edge (%d, %d)\n", rand_edge.first, rand_edge.second);
        tm.start();
        graph.insert_edge_update_muc(rand_edge.first, rand_edge.second);
        tm.stop();
        edge_stat.edge_ins_time = tm.interval();

        muc_id_t muc_id = graph.node_to_muc_vec[rand_edge.first];
        edge_stat.muc_num_edges = graph.muc_vec[muc_id].muc_subgraph.edge_set.size();
        edge_stat.muc_num_nodes = graph.muc_vec[muc_id].muc_subgraph.nodes_map.size();
        
        //estimating ideal speedup using the node and edge proportions
        double ideal_speedup;
        double V  = graph.size(); double E = graph.edge_set.size();
        double Vm = edge_stat.muc_num_nodes; double Em = edge_stat.muc_num_edges;
        ideal_speedup = ((V+E)/(Vm+Em))*(V/Vm);
        tot_ideal_speedup += ideal_speedup;
        ////
        
        printf("\tMUC # nodes: \t\t%d (%f G_n)\n",
                edge_stat.muc_num_nodes,
                edge_stat.muc_num_nodes/(double)graph.size());
        printf("\tMUC # edges: \t\t%d (%f G_M)\n",
                edge_stat.muc_num_edges,
                edge_stat.muc_num_edges/(double)graph.edge_set.size());
        printf("\n");
        printf("\tSpeedup (ideal): \t%f\n", ideal_speedup);
        printf("\n");
        printf("\tTime to update MUCs: \t%f\n", edge_stat.edge_ins_time);

        edge_stat_vec.push_back(edge_stat);
    }
    
    double avg_ideal_speedup = tot_ideal_speedup/num_edges_to_insert;
    printf("[[Average speedup (ideal): %f]]\n",
        avg_ideal_speedup);
}

/*
 * finds the fraction of nodes and edges in the largest bcc,
 * largest bcc in terms of number of nodes
 */
void lbcc_stat(graph_t& graph)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    vector<graph_hash_t> bcc_vec;
    graph.find_bicon_comp(bcc_vec);
    
    int lbcc_i = 0;//index of the lbcc, # of nodes
    int lbcc_j = 0;//index of the lbcc, # of edges
    int lbcc_n = 0;//highest # of nodes
    int lbcc_m = 0;//highest # of edges
    if(bcc_vec.size() > 1) {
        for(int i = 0; i < bcc_vec.size(); ++i) {
            if(bcc_vec[i].size() > lbcc_n) {
                lbcc_n = bcc_vec[i].size();
                lbcc_i = i;
            }
            if(bcc_vec[i].edge_set.size() > lbcc_m) {
                lbcc_m = bcc_vec[i].edge_set.size();
                lbcc_j = i;
            }
        }
    }
    if(lbcc_i != lbcc_j) {
        printf("LBCC_m is not LBCC_n!\n");
    }
    int num_nodes = lbcc_n;
    int num_edges = bcc_vec[lbcc_i].edge_set.size();
    double frac_nodes = (double)num_nodes/graph.size();
    double frac_edges = (double)num_edges/graph.edge_set.size();
    printf("\tLBCC # nodes:\t\t%d\n", num_nodes);
    printf("\tLBCC # edges:\t\t%d\n", num_edges);
    printf("\tLBCC frac nodes:\t%f\n", frac_nodes);
    printf("\tLBCC frac edges:\t%f\n", frac_edges);
    
}
