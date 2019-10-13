//
// Created by fanarosss on 12/10/19.
//
#include "compute_G.h"

using namespace std;

string amplify(vector<int>* H, int k){
    //amplifies h for k
    //g(x) = [h1(x)|h2(x)|h3(x)....|hk(x)].
    // TODO: see malloc sizes!

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    uniform_int_distribution<int> distribution (0, H->size());
    string gx;

    for (int j = 0; j < k; j++) {
        gx.append(to_string((*H)[distribution(generator)]));
    }

    return gx;
}
