

#include "utility.h"
#include <algorithm>
#include <numeric>
#include <math.h>
#include <iostream>

using namespace std;

void simple_stats(vector<double> vec,
        double& mean,
        double& stddev,
        double& median,
        double& min,
        double& max)
{
    sort(vec.begin(), vec.end());
    min = *(min_element(vec.begin(), vec.end()));
    max = *(max_element(vec.begin(), vec.end()));
    double sum = accumulate(vec.begin(), vec.end(), 0.0);
    mean = sum/vec.size();
    median = vec[vec.size()/2];
    
    double sum_2 = 0;
    for(int i = 0; i < vec.size(); ++i)
        sum_2 += ((vec[i]-mean)*(vec[i]-mean));
    stddev = sqrt(sum_2/vec.size());
    
}


string extract_graph_name(string path)
{
    string out;
    int pos = path.rfind("/") + 1;
    int length = path.length() - pos + 1;
    out = path.substr(pos, length);
    return out;
}