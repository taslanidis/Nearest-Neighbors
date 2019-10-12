//
// Created by fanarosss on 12/10/19.
//
#include "compute_G.h"

using namespace std;

string amplify(vector<int>* H){
    //amplifies h for k
    //g(x) = [h1(x)|h2(x)|h3(x)....|hk(x)].
    // TODO: see malloc sizes!

    string gx;
    for (int i = 0; i < H->size(); i++) {
        gx.append(to_string((*H)[i]));
    }

    return gx;
}
