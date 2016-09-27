#include <algorithm>
#include <cmath>

#include "experiments.h"
#include "stdlib.h"
#include "utility.h"
#include "bc.h"


double  exp_brandes_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double time_to_run_all_iter)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    
    double  brandes_time;
    bool    do_all_iter = false;
    timer   tm;
    
    printf("Brandes runtime: \t");
    if(num_iter != -1) {
        tm.start();
        tr1_map_t(double) bc_vec = brandes_bc_hash_out(graph, num_iter);
        tm.stop();
        brandes_time = (tm.interval()*graph.size()/num_iter);
        if(brandes_time < time_to_run_all_iter) {
            do_all_iter = true;
        }
    } else {
        do_all_iter = true;
    }

    if(do_all_iter) {
        num_iter = -1;
        tm.start();
        tr1_map_t(double) bc_vec = brandes_bc_hash_out(graph, num_iter);
        tm.stop();
        brandes_time = tm.interval();
        printf("%f\n", brandes_time);
        printf("Brandes avg iter time: \t%f\n", tm.interval()/graph.size());
        printf("DID ALL ITERATIONS...\n");
    } else {
        printf("%f\n", brandes_time);
        printf("Brandes avg iter time: \t%f\n", tm.interval()/num_iter);
        printf("DID %d ITERATIONS...\n", num_iter);
    }
    
    return brandes_time;
}


void    exp_inc_brandes_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double brandes_time)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    timer tm;
    vector<edge_stat_t> edge_stat_vec;
    
    vector<double> bc_vec;
    bc_vec.resize(graph.size());
 
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
    
    for(int i = 0; i < cnt_vec2.size(); ++i) {
        tot_cnt_arr[0] += cnt_vec2[i][0];
        tot_cnt_arr[1] += cnt_vec2[i][1];
        tot_cnt_arr[2] += cnt_vec2[i][2];
        
        tot_time_arr[0] += time_vec2[i][0];
        tot_time_arr[1] += time_vec2[i][1];
        tot_time_arr[2] += time_vec2[i][2];
    }
    double tot_speedup = 0;
    for(int i = 0; i < cnt_vec2.size(); ++i) {
        double inc_br_time = time_vec2[i][0] + time_vec2[i][1] + time_vec2[i][2];
        tot_speedup += (brandes_time/inc_br_time);
    }
    
    double inc_tot_time = tot_time_arr[0] + tot_time_arr[1] + tot_time_arr[2];
    double inc_avg_time = inc_tot_time/rand_edge_vec.size();
    printf("Incremental-Brandes avg. runtime:\t%f\n", inc_avg_time);
    double avg_d1_tm = tot_time_arr[1]/tot_cnt_arr[1];
    double avg_d2_tm = tot_time_arr[2]/tot_cnt_arr[2];
    double avg_br_tm = brandes_time/graph.size();
    //the prints below is for the total iterations for all the inserted edges
    //they don't make sense if not reported per edge
//    printf("Tot d0 BFS:\t\t\t\t%d/%d\n", tot_cnt_arr[0], graph.size());
//    printf("Tot d1 BFS:\t\t\t\t%d/%d\n", tot_cnt_arr[1], graph.size());
//    printf("Tot d2 BFS:\t\t\t\t%d/%d\n", tot_cnt_arr[2], graph.size());
    printf("Avg d1 BFS time:\t\t\t%f (%f)\n", avg_d1_tm, avg_d1_tm/avg_br_tm);
    printf("Avg d2 BFS time:\t\t\t%f (%f)\n", avg_d2_tm, avg_d2_tm/avg_br_tm);
    printf("[[Avg. speedup (inc_brandes): %f]]\n", brandes_time/inc_avg_time);
    printf("[[Avg. speedup (inc_brandes)(proper): %f]]\n", tot_speedup/rand_edge_vec.size());
}


void    exp_qube_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double brandes_time)
{
//    printf("==========================================================\n");
//    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
//    printf("Graph # nodes: \t\t%d\n", graph.size());
//    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    printf("Graph[%s]  V[%d]  E[%d]  Brandes_tm[%.2f]\n",
        graph.graph_name.c_str(),
        graph.size(),
        graph.edge_set.size(),
        brandes_time);
    
    double muc_find_time;
    timer tm;
    tr1_map_t(double) bc_vec = brandes_bc_hash_out(graph, 1);
    printf("# MUCs: \t\t");
    tm.start();
    graph.find_mucs();
    tm.stop();
    muc_find_time = tm.interval();
    printf("%d\n", graph.get_num_mucs());
    printf("Time to find MUCs: \t%f\n", muc_find_time);
    
//    if(num_iter == -1) {
//        printf("WILL DO ALL ITERATIONS...\n");
//    } else {
//        printf("WILL DO %d ITERATIONS...\n", num_iter);
//    }
    
    vector<edge_stat_t> edge_stat_vec;
    
    double tot_speedup = 0;
    double tot_time = 0;
    
    vector<double> speedup_vec;
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        edge_t rand_edge = rand_edge_vec[i];
        edge_stat_t edge_stat;
//        printf("Adding edge (%d, %d)\n", rand_edge.first, rand_edge.second);
        tm.start();
        graph.insert_edge_update_muc(rand_edge.first, rand_edge.second);
        tm.stop();
        edge_stat.edge_ins_time = tm.interval();

        muc_id_t muc_id = graph.node_to_muc_vec[rand_edge.first];
        edge_stat.muc_num_edges = graph.muc_vec[muc_id].muc_subgraph.edge_set.size();
        edge_stat.muc_num_nodes = graph.muc_vec[muc_id].muc_subgraph.nodes_map.size();
        tr1::unordered_map<node_id_t, double> bc_map;
        tm.start();
        graph.muc_vec[muc_id].compute_bc(bc_map, num_iter);
        tm.stop();
        if(num_iter == -1)
            edge_stat.muc_bc_update_time = tm.interval();
        else
            edge_stat.muc_bc_update_time = (tm.interval()/num_iter)*edge_stat.muc_num_nodes;
        //TMP CODE XXX FIX THIS SHIT
        edge_stat.muc_bc_update_time += edge_stat.edge_ins_time;
        ////
        
        double speedup = brandes_time/edge_stat.muc_bc_update_time;
        tot_speedup += speedup;
        speedup_vec.push_back(speedup);
        
        tot_time += edge_stat.muc_bc_update_time;
        
//        printf("\tMUC # nodes: \t\t%d (%f G_n)\n",
//                edge_stat.muc_num_nodes,
//                edge_stat.muc_num_nodes/(double)graph.size());
//        printf("\tMUC # edges: \t\t%d (%f G_M)\n",
//                edge_stat.muc_num_edges,
//                edge_stat.muc_num_edges/(double)graph.edge_set.size());
//        printf("\n");
//        printf("\tTime to update BC: \t%f\n", edge_stat.muc_bc_update_time);
//        printf("\tSpeedup: \t\t%f\n", speedup);
//        printf("\n");
//        printf("\tQUBE avg iter time: \t%f\n",
//                edge_stat.muc_bc_update_time/edge_stat.muc_num_nodes);
//        printf("\tTime to update MUCs: \t%f\n", edge_stat.edge_ins_time);
        
        printf("e(%-6d,%-6d)  tm[%.2f]  sup[%.2f]\n",
                rand_edge.first,
                rand_edge.second,
                edge_stat.muc_bc_update_time,
                speedup);

        edge_stat_vec.push_back(edge_stat);
    }
    
    double avg_speedup = tot_speedup/rand_edge_vec.size();
//    printf("[[Avg. speedup (qube): %f]]\n", avg_speedup);
    double avg_time = tot_time/rand_edge_vec.size();
    printf("Avg.tm[%.2f]  Avg.sup[%.2f]\n\n", avg_time, avg_speedup);
}


void    exp_inc_qube_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double brandes_time)
{
//    printf("==========================================================\n");
//    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
//    printf("Graph # nodes: \t\t%d\n", graph.size());
//    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    printf("Graph[%s]  V[%d]  E[%d]  Brandes_tm[%.2f]\n",
        graph.graph_name.c_str(),
        graph.size(),
        graph.edge_set.size(),
        brandes_time);
    double muc_find_time;
    timer tm;
    
    tr1_map_t(double) bc_vec = brandes_bc_hash_out(graph, 1);
    
    printf("# MUCs: \t\t");
    tm.start();
    graph.find_mucs();
    tm.stop();
    muc_find_time = tm.interval();
    printf("%d\n", graph.get_num_mucs());
    printf("Time to find MUCs: \t%f\n", muc_find_time);
    
//    if(num_iter == -1) {
//        printf("WILL DO ALL ITERATIONS...\n");
//    } else {
//        printf("WILL DO %d ITERATIONS...\n", num_iter);
//    }
    
    vector<edge_stat_t> edge_stat_vec;
    
    double tot_speedup = 0;
    double tot_time = 0;
    
    vector<double> speedup_vec;
    for(int i = 0; i < rand_edge_vec.size(); ++i) {
        edge_t rand_edge = rand_edge_vec[i];
        edge_stat_t edge_stat;
//        printf("Adding edge (%d, %d)\n", rand_edge.first, rand_edge.second);
        tm.start();
        graph.insert_edge_update_muc(rand_edge.first, rand_edge.second);
        tm.stop();
        edge_stat.edge_ins_time = tm.interval();

        muc_id_t muc_id = graph.node_to_muc_vec[rand_edge.first];
        edge_stat.muc_num_edges = graph.muc_vec[muc_id].muc_subgraph.edge_set.size();
        edge_stat.muc_num_nodes = graph.muc_vec[muc_id].muc_subgraph.nodes_map.size();
        tr1::unordered_map<node_id_t, double> bc_map;
        tm.start();
        graph.muc_vec[muc_id].compute_bc_inc(bc_map, rand_edge.first, rand_edge.second, num_iter);
        tm.stop();
        if(num_iter == -1)
            edge_stat.muc_bc_update_time = tm.interval();
        else
            edge_stat.muc_bc_update_time = (tm.interval()/num_iter)*edge_stat.muc_num_nodes;;
        //TMP CODE XXX FIX THIS SHIT
        edge_stat.muc_bc_update_time += edge_stat.edge_ins_time;
        ////
        
        double speedup = brandes_time/edge_stat.muc_bc_update_time;
        tot_speedup += speedup;
        speedup_vec.push_back(speedup);
        
        tot_time += edge_stat.muc_bc_update_time;
        
//        printf("\tMUC # nodes: \t\t%d (%f G_n)\n",
//                edge_stat.muc_num_nodes,
//                edge_stat.muc_num_nodes/(double)graph.size());
//        printf("\tMUC # edges: \t\t%d (%f G_M)\n",
//                edge_stat.muc_num_edges,
//                edge_stat.muc_num_edges/(double)graph.edge_set.size());
//        printf("\n");
//        printf("\tTime to update BC: \t%f\n", edge_stat.muc_bc_update_time);
//        printf("\tSpeedup: \t\t%f\n", speedup);
//        printf("\n");
//        printf("\tQUBE avg iter time: \t%f\n",
//                edge_stat.muc_bc_update_time/edge_stat.muc_num_nodes);
//        printf("\tTime to update MUCs: \t%f\n", edge_stat.edge_ins_time);
        printf("e(%-6d,%-6d)  tm[%.2f]  sup[%.2f]\n",
                rand_edge.first,
                rand_edge.second,
                edge_stat.muc_bc_update_time,
                speedup);

        edge_stat_vec.push_back(edge_stat);
    }
    
    double avg_speedup = tot_speedup/rand_edge_vec.size();
//    printf("[[Avg. speedup (inc_qube): %f]]\n", avg_speedup);
    double avg_time = tot_time/rand_edge_vec.size();
    printf("Avg.tm[%.2f]  Avg.sup[%.2f]\n\n", avg_time, avg_speedup);
}


void    exp_fuad_p(graph_t& graph,
            int num_iter,
            vector<edge_t> rand_edge_vec,
            double brandes_time)
{
    printf("==========================================================\n");
    printf("Graph: \t\t\t%s\n", graph.graph_name.c_str());
    printf("Graph # nodes: \t\t%d\n", graph.size());
    printf("Graph # edges: \t\t%d\n", graph.edge_set.size());
    
    tr1_map_t(double) bc_vec = brandes_bc_hash_out(graph, 1);
    timer tm;

    double tot_fuad_time = 0;
    double tot_speedup = 0;
    if(num_iter == -1) {
        printf("WILL DO ALL ITERATIONS...\n");
    } else {
        printf("WILL DO %d ITERATIONS...\n", num_iter);
    }
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
        if(num_iter == -1) {
            est_time = bcc_stat.bc_update_time + bcc_stat.bcc_find_time;
        }
        
        double avg_iter_time_naive;
        graph.approx_bcc_iter_tm(rand_edge.first, rand_edge.second, avg_iter_time_naive, 10);
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
        if(num_iter != -1) {
            est_time += bcc_stat.tot_d1_tm*bcc_stat.tot_d1_iter/bcc_stat.num_d1_iter;
            est_time += bcc_stat.tot_d2_tm*bcc_stat.tot_d2_iter/bcc_stat.num_d2_iter;
        }
       
        tot_speedup += (brandes_time/est_time);
        speedup_vec.push_back(brandes_time/est_time);
        
        printf("\tBCC num nodes: \t\t\t%d (%f G_n)\n",
                bcc_stat.bcc_num_nodes,
                (double)bcc_stat.bcc_num_nodes/graph.size());
        printf("\tBCC num edges: \t\t\t%d (%f G_m)\n",
                bcc_stat.bcc_num_edges,
                (double)bcc_stat.bcc_num_edges/graph.edge_set.size());
        printf("\n");
        printf("\tTime to update BC: \t\t%f\n",
                est_time);
        printf("\tSpeedup: \t\t\t%f\n",
                brandes_time/est_time);
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
    double avg_speedup = tot_speedup/rand_edge_vec.size();
    
    printf("[[Avg. speedup (fuad): %f]]\n", avg_speedup);
    
}