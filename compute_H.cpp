//
// Created by kostas on 12/10/19.
//

#include "compute_H.h"

using namespace std;

void compute_H (vector<int>* H, vector<vector<int>> a, int d, int k){
    int m, M;
    double power = pow(2,32);
    m = (int)(power - 5);
    M = pow(2, 32/k);
    cout << "m is: " << m << endl;
    cout << "M is: " << M << endl;
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
