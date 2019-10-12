//
// Created by kostas on 12/10/19.
//

#include "compute_s.h"

using namespace std;

void compute_s(vector<double>* s, int w, int d){
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    uniform_real_distribution<double> distribution (0.0, w);

    for (int i=0; i<d; ++i) {
        s->push_back(distribution(generator));
    }
}