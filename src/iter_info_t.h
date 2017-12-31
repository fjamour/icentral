/* 
 * File:   iter_info_t.h
 * Author: fuad
 *
 * Created on November 24, 2014, 4:09 PM
 */


#ifndef ITER_INFO_T_H
#define	ITER_INFO_T_H

#include <vector>
#include "types.h"
using namespace std;

/*
 * Stores the information used in one iteration of iCentral
 * Namely, \sigma \P \delta \D \S
 * both before edge insertion and after edge insertion
 * TODO study the effectiveness of storing old/new
 * TODO reconsider what to store and what not to
 *      (after the iCentral is implemented)
 * TODO reconsider the names of the structures
 */
struct iter_info_t {
    vector<vector<node_id_t> >       P;
    vector<int>                      sigma_vec;
    vector<int>                      dist_vec;
    vector<double>                   delta_vec;
    vector<double>                   Delta_vec;
    
    vector<int>                      new_sigma_vec;
    vector<double>                   new_delta_vec;
    vector<double>                   new_Delta_vec;

    vector<int>                      sigma_inc_vec;
    
    vector<bool>                     visited_vec;
    vector<node_id_t>                S;
    
    void init_all(size_t N);
    void init_new(size_t N);
};


#endif	/* ITER_INFO_T_H */

