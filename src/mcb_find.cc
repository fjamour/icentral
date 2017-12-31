#include <iostream>

#include "mcb_find.h"

#ifndef NO_LEDA
#include <LEP/mcb/mcb.h>
#endif

#include "graph_t.h"

using namespace std;

#ifndef NO_LEDA
void fill_leda_graph(graph_t* graph, leda::graph* leda_g )
{
    vector<leda::graph::node> node_vec;
    for(int i = 0; i < graph->size(); i++) {
        node_vec.push_back(leda_g->new_node());
    }
    for(set<pair<node_id_t, node_id_t> >::iterator it = graph->edge_set.begin(); it != graph->edge_set.end(); it++) {
        leda_g->new_edge(node_vec[it->first], node_vec[it->second]);
    }
    leda_g->make_undirected();
}
#endif

void mcb_find(graph_t* graph, mcb_t* mcb_out)
{
#ifndef NO_LEDA
    leda::graph G;
    fill_leda_graph(graph, &G);

    leda::edge_array<int> len(G, 1);
    mcb::edge_num enumb(G);
    leda::array< mcb::spvecgf2 > mcb;
    int weight = mcb::UMCB_SVA(G, len, mcb, enumb);

    //cout << "Number of nodes: " << G.number_of_nodes() << endl;
    //cout << "Number of edges: " << G.number_of_edges() << endl;
    //G.print(cout);

    int i, j;
    leda::edge e;
    for (i = 0; i < enumb.dim_cycle_space(); ++i) {
        //printf("Cycle: %d\n", i);
        cycle_t cycle;
        forall(j, mcb[i]) { // traverse edges of i-th cycle
            e = enumb(j);
            // do something with edge e
            //G.print_edge(e);
            //cout << endl;
            node_id_t src = G.source(e)->id();
            node_id_t dst = G.target(e)->id();
            edge_t edge(src, dst);
            cycle.push_back(edge);
        }
        mcb_out->cycle_vec.push_back(cycle);
    }
#endif
}
