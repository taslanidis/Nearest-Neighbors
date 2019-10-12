//
// Created by kostas on 12/10/19.
//

#include "compute_H.h"

using namespace std;

void compute_H (vector<int>* H, vector<vector<int>> a, int d, int k, int w){
    int m, M;
    m = w;
    M = pow(2, 32/k);
    int h, term;
    for (int i=0; i<a.size(); i++){
        h=0;
        for (int j=0; j<d; j++){
            term = a[i][d-1-j]*pow(m,j);
            h += term % M;
        }
        H->push_back(h);
    }
}
