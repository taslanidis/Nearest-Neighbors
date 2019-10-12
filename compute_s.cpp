//
// Created by kostas on 12/10/19.
//

#include "compute_s.h"

using namespace std;

void compute_s(vector<int>* s, int w, int d){
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    uniform_int_distribution<int> distribution (0, w);

    for (int i=0; i<d; ++i) {
        s->push_back(distribution(generator));
    }
}