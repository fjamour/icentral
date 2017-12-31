/* 
 * File:   component_t.h
 * Author: fuad
 *
 * Created on November 24, 2014, 4:08 PM
 */

#ifndef COMPONENT_T_H
#define	COMPONENT_T_H

#include "graph_t.h"

enum comp_type_t {BCC, MUC, GRAPH};

/*
 * The component (subgraph along with other needed information)
 * that bc computation blocks operate on.
 * Could be a BCC, an MUC, or just a graph.
 */
struct component_t {
    //maps articulation points to sizes of subgraphs connected to the bcc
    //through them
    typedef tr1::unordered_map<node_id_t, vector<int> >   art_pt_map_t;
    
    subgraph_t      subgraph;
    art_pt_map_t    art_pt_map;
    comp_type_t     comp_type;
    
    void print()
    {
        printf("\n");
        subgraph.print_graph();
        art_pt_map_t::iterator it = art_pt_map.begin();
        for(; it != art_pt_map.end(); ++it) {
            printf("Articulation point [%d]\n", it->first);
            for(int i = 0; i < it->second.size(); ++i) {
                printf("\t%d", it->second[i]);
            }
            printf("\n");
        }
    }
};

#endif	/* COMPONENT_T_H */

