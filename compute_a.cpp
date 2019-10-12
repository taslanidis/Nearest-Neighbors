//
// Created by fanarosss on 12/10/19.
//
#include "compute_a.h"

using namespace std;

void compute_a(vector<vector<int>>* a_projects, vector<vector<int>> x, vector<int>* s, int w, int d) {
    // Ai = (Xi - Si) / W
    // Project every X to A in d-dimensional grid shifted by S, where every cell size = W

    int ai;
    vector<int> a;
    for (int i = 0; i < x.size(); i++) {
        for (int dim = 0; dim < d; dim++) {
            ai = floor((double)(x[i][dim] - (*s)[dim]) / w);
            a.push_back(ai);
        }
        a_projects->push_back(a);
        a.clear();
    }

}