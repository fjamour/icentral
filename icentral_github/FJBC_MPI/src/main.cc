#include <iostream>
#include <stdio.h>
#include <fstream>
#include <list>
#include <algorithm>
#include <cmath>


#include "graph_t.h"
#include "bc.h"
#include "experiments.h"
#include "utility.h"

#include "unit_tests.h"



using namespace std;


void do_paper_exp(int num_edges,
        int num_iter,
        double max_time,
        int rand_seed,
        vector<string> path_vec,
        bool do_inc_brandes,
        bool do_qube,
        bool do_inc_qube,
        bool do_fuad) 
{     
    //to make sure the exact same edges are used for my system and QUBE
    //I pass the same edges to both systems
    
    vector<double>  brandes_time_vec;
    
    vector<vector<edge_t> > edge_vec2;
    printf("Reading graphs and generating edges...\n");
    printf("========================================\n");
    for(int i = 0; i < path_vec.size(); ++i) {
        graph_t graph;
        string path = path_vec[i].c_str();
        graph.read_graph(path);
        graph.graph_name = extract_graph_name(path_vec[i]);
        vector<edge_t> edge_vec;
        gen_rand_edges(num_edges, graph, edge_vec);
        edge_vec2.push_back(edge_vec);
    }
    
    printf("\n\n\n");
    printf("Starting Brandes...\n");
    printf("========================================\n");
    for(int i = 0; i < path_vec.size(); ++i) {
        graph_t graph;
        string path = path_vec[i].c_str();
        graph.read_graph(path);
        graph.graph_name = extract_graph_name(path_vec[i]);
        double brandes_time;
        brandes_time = exp_brandes_p(graph, num_iter, edge_vec2[i], max_time);
        brandes_time_vec.push_back(brandes_time);
    }
    
    if(do_inc_brandes) {
        printf("\n\n\n");
        printf("Starting Incremental-Brandes...\n");
        printf("========================================\n");
        for(int i = 0; i < path_vec.size(); ++i) {
            graph_t graph;
            string path = path_vec[i].c_str();
            graph.read_graph(path);
            graph.graph_name = extract_graph_name(path_vec[i]);
            double brandes_time = brandes_time_vec[i];
            exp_inc_brandes_p(graph, num_iter, edge_vec2[i], brandes_time);
        }
    }
    
    if(do_qube) {
        printf("\n\n\n");
        printf("Starting QUBE...\n");
        printf("========================================\n");
        for(int i = 0; i < path_vec.size(); ++i) {
            graph_t graph;
            string path = path_vec[i].c_str();
            graph.read_graph(path);
            graph.graph_name = extract_graph_name(path_vec[i]);
            double brandes_time = brandes_time_vec[i];
            int loc_num_iter = num_iter;
            if(num_iter != -1 && brandes_time < max_time)
                loc_num_iter = -1;
            exp_qube_p(graph, loc_num_iter, edge_vec2[i], brandes_time);
        }
    }
    
    if(do_inc_qube) {
        printf("\n\n\n");
        printf("Starting Incremental-QUBE...\n");
        printf("========================================\n");
        for(int i = 0; i < path_vec.size(); ++i) {
            graph_t graph;
            string path = path_vec[i].c_str();
            graph.read_graph(path);
            graph.graph_name = extract_graph_name(path_vec[i]);
            double brandes_time = brandes_time_vec[i];
            int loc_num_iter = num_iter;
            if(num_iter != -1 && brandes_time < max_time)
                loc_num_iter = -1;
            exp_inc_qube_p(graph, loc_num_iter, edge_vec2[i], brandes_time);
        }
    }
    
    if(do_fuad) {
        printf("\n\n\n");
        printf("Starting FUAD...\n");
        printf("========================================\n");
        for(int i = 0; i < path_vec.size(); ++i) {
            graph_t graph;
            string path = path_vec[i].c_str();
            graph.read_graph(path);
            graph.graph_name = extract_graph_name(path_vec[i]);
            double brandes_time = brandes_time_vec[i];
            int loc_num_iter = num_iter;
            if(num_iter != -1 && brandes_time < max_time)
                loc_num_iter = -1;
            exp_fuad_p(graph, loc_num_iter, edge_vec2[i], brandes_time);
        }
    }
}

void paper_exp_main(int argc, char** argv)
{
    int             num_edges, num_iter, rand_seed;
    bool            do_inc_brandes, do_qube, do_inc_qube, do_fuad;
    double          max_time;
    vector<string>  path_vec;
    
    if(argc != 2) {
        printf("Pass one parameter, path with experiment details\n");
        printf("num_edges, num_iter, max_time, rand_seed\n");
        printf("do_inc_brandes, do_qube, do_inc_qube, do_fuad\n");
        printf("list of graph paths\n");
        exit(1);
    } else {
        FILE* fin = fopen(argv[1], "r");
        fscanf(fin, "%d, %d, %lf, %d\n", &num_edges, &num_iter, &max_time, &rand_seed);
        int t1, t2, t3, t4;
        fscanf(fin, "%d, %d, %d, %d\n", &t1, &t2, &t3, &t4);
        do_inc_brandes  = (t1 != 0);
        do_qube         = (t2 != 0);
        do_inc_qube     = (t3 != 0);
        do_fuad         = (t4 != 0);
        char buff[1024*4];
        while(fscanf(fin, "%s\n", buff) != EOF) {
            string path = buff;
            path_vec.push_back(path);
        }
    }
    do_paper_exp(num_edges, num_iter, max_time, rand_seed, path_vec,
        do_inc_brandes, do_qube, do_inc_qube, do_fuad);
}



vector<string> fill_path_vec()
{
    vector<string> vec;
    string path_f  = "/home/jamourft/Desktop/Research/Betweenness-Centrality/data/fj_lcc_graphs/";
    string graphs_arr[] = 
        {   
          "Erdos02.lcc.net"
        , "Erdos972.lcc.net"
        , "Cagr.lcc.net"
        , "Eva.lcc.net"
        , "Epa.lcc.net"
        , "Contact.lcc.net"
        , "Wiki-Vote.lcc.net"
        };
    
    for(int i = 0; i < 7; ++i) {
        string s = path_f + graphs_arr[i];
        vec.push_back(s);
    }
    return vec;
}



void kdd_exp_main(int argc, char** argv, int rank, int size)
{
    int                     num_edges, rand_seed, num_threads;
    bool                    do_icent, do_bcc_icent, ext_edges;
    bool                    do_fast_brandes, do_brandes, do_qube, do_inc_qube;
    operation_t             op;
    //for qube and incremental qube, old implementation will be used, and
    //no parallelization
    vector<string>          path_vec;
    vector<vector<edge_t> > edge_vec2;
    vector<double>          brandes_tm_vec;
    
    MPI_Status    status;
    
    if(argc != 2) {
        if(rank == 0) {
            printf("Pass one parameter, path with experiment details\n");
            printf("deletion\n");
            printf("num_edges, num_threads, rand_seed\n");
            printf("external_edges, do_icent, do_bcc_icent\n");
            printf("do_fast_brandes, do_brandes, do_qube, do_inc_qube\n");
            printf("list of graph paths\n");
            printf("if external_edges is nonzero a file with graph_name.edges is expected\n");
            printf("external edges are assumed to be in the graph\n");
            printf("external edges should not be bridges\n");
            printf("if deletion is 1, bcc_icent with deletions will be invoked");
        }
        return;
        //exit(1);
    } else {
        FILE* fin = fopen(argv[1], "r");
        int del_int;
        fscanf(fin, "%d", &del_int);
        if(del_int == 1) {
            op = DELETION;
        } else {
            op = INSERTION;
        }
        fscanf(fin, "%d, %d, %d\n", &num_edges, &num_threads, &rand_seed);
        int t1, t2, t3, t4;
        fscanf(fin, "%d, %d, %d\n", &t1, &t2, &t3);
        ext_edges       = (t1 != 0);
        do_icent        = (t2 != 0);
        do_bcc_icent    = (t3 != 0);
        fscanf(fin, "%d, %d, %d, %d\n", &t1, &t2, &t3, &t4);
        do_fast_brandes = (t1 != 0);
        do_brandes      = (t2 != 0);
        do_qube         = (t3 != 0);
        do_inc_qube     = (t4 != 0);
        
        char buff[1024*4];
        while(fscanf(fin, "%s\n", buff) != EOF) {
            string path = buff;
            path_vec.push_back(path);
        }
        fclose(fin);
        //fill edges, read from file if needed
        if(ext_edges) {
            //these are edges in the graph, and have to be removed then inserted
            //in order
            for(int i = 0; i < path_vec.size(); ++i) {
                vector<edge_t> edge_vec;
                string edge_file_path = path_vec[i] + ".edges";
                fin = fopen(edge_file_path.c_str(), "r");
                int v1, v2;
                while(fscanf(fin, "%d %d\n", &v1, &v2) != EOF) {
                    edge_vec.push_back(make_pair(v1, v2));
                }
                edge_vec2.push_back(edge_vec);
                fclose(fin);
            }
        } else {
            for(int p = 0; p < path_vec.size(); ++p) {
                string graph_path = path_vec[p];
                graph_t graph;
                graph.read_graph(graph_path);
                graph.graph_name = extract_graph_name(graph_path);
                vector<edge_t> edge_vec;
                
                //gen_rand_edges(num_edges, graph, edge_vec);
                //edge_vec2.push_back(edge_vec);
                if(rank == 0) {
                    srand(rand_seed);
                    //master generates random edges and sends to everyone
                    if(op == INSERTION) {
                        gen_rand_edges(num_edges, graph, edge_vec);
                    } else if(op == DELETION) {
                        gen_rand_edges_deletions(num_edges, graph, edge_vec);
                    }
                    
                    //TMP code to make sure new random edges are generated.. Safely remove!
                    //for(int re = 0; re < edge_vec.size(); ++re) {
                    //    printf("e(%d, %d)\n", edge_vec[re].first, edge_vec[re].second);
                    //}
                    /////
                    //printf("RANK[%d] -- e[0]: %u %u\n", rank, edge_vec[4].first, edge_vec[4].second);
                    for(int p = 1; p < size; ++p) {
                        MPI_Send(&edge_vec[0], edge_vec.size()*sizeof(edge_t), MPI_CHAR, p, 0, MPI_COMM_WORLD);
                    }
                } else {
                    //slaves get random edges from the master
                    edge_vec.resize(num_edges);
                    MPI_Recv(&edge_vec[0], edge_vec.size()*sizeof(edge_t), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    //printf("RANK[%d] -- e[0]: %u %u\n", rank, edge_vec[4].first, edge_vec[4].second);
                }
                edge_vec2.push_back(edge_vec);
            }
        }
    }
    
    int num_graphs = path_vec.size();
    brandes_tm_vec.resize(num_graphs);
    fill(brandes_tm_vec.begin(), brandes_tm_vec.end(), 1.0);
    
    timer tm;
    vector<double> BC_vec;
    double brandes_time;
    
    if(do_bcc_icent) {
        if(rank == 0) {
            printf("\n\n\n");
            printf("Starting BCC+iCentral [%d threads] [%s]...\n", num_threads, (op==DELETION?"DELETION":"INSERTION"));
            printf("========================================\n");
        }
        for(int i = 0; i < path_vec.size(); ++i) {
            graph_t graph;
            string path = path_vec[i].c_str();
            graph.read_graph(path);
            graph.graph_name = extract_graph_name(path_vec[i]);
            timing__Update_BC_graph(
                graph,
                edge_vec2[i], 
                -1,//do all sources (not approximation) 
                BCC, 
                false, 
                brandes_tm_vec[i], 
                num_threads,
                ext_edges,
                op
                );
            //synchronization barrier so that no one starts next graph before others
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    

    

}


//int main(int argc, char**argv) {
//    setbuf(stdout, NULL);
//    cout << "Hello!" << endl;
//
//    
//    //paper_exp_main(argc, argv);
//    
//    //test__brandes_BC();
//    //test__Update_BC();
//    //test__iCentral();
//    
////    vector<string> path_vec = fill_path_vec();
////    timing__Update_BC(path_vec, 10, 111, -1, GRAPH, 4);
//    
////    vector<string> path_vec = fill_path_vec();
////    timing__Update_BC_mem(path_vec, 10, 111, -1, 111);
//    
////    test__Update_BC_mem();
//    
//    
//    kdd_exp_main(argc, argv);
//    return 0;
//}

#include <stdio.h>
#include <string.h>
#include <mpi.h>


/*
 * 
 */
void run_parallel_brandes(
        string  graph_path,
        string  output_paht,
        int     num_threads
        )
{
   //1. read input graph
   //2. call run parallel brandes
   //3. write the result to disk (label: BC value)
}

int main( int argc, char *argv[] )
{
    int i;
    int rank;
    int size;
    MPI_Status    status;
    char str_message[100];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

//    sprintf(str_message, "Hello world from process %d of %d\n", rank, size);
//    int cnt=strlen(str_message)+1;
//    if(rank!=0) MPI_Send(str_message, cnt, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
//    if(rank==0)
//    {
//        for(i=1;i<size;i++)
//        {
//            MPI_Recv(str_message, 100, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//
//            // If we need to receive messages in order, we should do this instead:
//            //MPI_Recv(str_message, 100, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//
//            printf("%s", str_message);
//        }
//    }

//    if(rank == 0) {
//        vector<double> vv;
//        vv.push_back(10);
//        vv.push_back(20);
//        MPI_Send(&vv[0], vv.size(), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
//        
//        vector<pair<unsigned int, unsigned int> > vvp;
//        vvp.push_back(make_pair(20, 17));
//        vvp.push_back(make_pair(25, 27));
//        MPI_Send(&vvp[0], vvp.size()*sizeof(pair<unsigned int, unsigned int>), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
//    }
//    if(rank == 1) {
//        vector<double> vv;
//        vv.resize(2);
//        MPI_Recv(&vv[0], vv.size(), MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//        printf("%f, %f\n", vv[0], vv[1]);
//        vector<pair<unsigned int, unsigned int> > vvp;
//        vvp.resize(2);
//        MPI_Recv(&vvp[0], vvp.size()*sizeof(pair<unsigned int, unsigned int>), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//        printf("%u %u %u %u\n", vvp[0].first, vvp[0].second, vvp[1].first, vvp[1].second);
//        
//    }
    
    
    kdd_exp_main(argc, argv, rank, size);

    MPI_Finalize();
    return 0;
}