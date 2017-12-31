#include "iter_info_t.h"
#include "utility.h"

void iter_info_t::init_all(size_t N)
{
    fill_vec<vector<node_id_t> >(P, N, vector<node_id_t>());
    fill_vec<int>(sigma_vec, N, 0);
    fill_vec<int>(dist_vec, N, -1);
    fill_vec<int>(sigma_inc_vec, N, 0);
    
    fill_vec<double>(delta_vec, N, 0);
    fill_vec<double>(Delta_vec, N, 0);

    fill_vec<bool>(visited_vec, N, false);
    S.clear();
    
//    fill_vec<int>(new_sigma_vec, N, 0);
//    fill_vec<double>(new_delta_vec, N, 0);
//    fill_vec<double>(new_Delta_vec, N, 0);
    
//    fill_vec<vector<node_id_t> >(old_P, N, vector<node_id_t>());
//    fill_vec<int>(old_sigma_vec, N, 0);
//    fill_vec<int>(old_dist_vec, N, -1);
//    fill_vec<double>(old_delta_vec, N, 0);
//    fill_vec<double>(old_Delta_vec, N, 0);
}

//needed for d=1 only
void iter_info_t::init_new(size_t N)
{
    fill_vec<int>(new_sigma_vec, N, 0);
    fill_vec<double>(new_delta_vec, N, 0);
    fill_vec<double>(new_Delta_vec, N, 0);
}