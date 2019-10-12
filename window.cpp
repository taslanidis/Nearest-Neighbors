//
// Created by fanarosss on 12/10/19.
//
#include "window.h"

using namespace std;

int compute_window(vector<vector<int>> dataset) {
    // 1. take all points in dataset
    // 2. find their nearest neighbor using L1 metric
    // 3. average all distances between points
    // L1 = sum(|P1i - P2i|) for i in 0,d-1

    vector<int> distances;
    int L1, min_distance;
    int size = dataset.size();
    int d = dataset[0].size();
    vector<int> P1;
    vector<int> P2;
    for (int i = 0; i < size; i++) {
        P1 = dataset[i];
        L1 = 0;
        min_distance = -1;
        for (int j = 0; j < size; j++) {
            if (i != j) {
                P2 = dataset[j];
                for (int dim = 0; dim < d; dim++)
                    L1 += abs(P1[dim] - P2[dim]);
                if (min_distance == -1)
                    min_distance = L1;
                if (L1 < min_distance)
                    min_distance = L1;
            }
            L1 = 0;
        }
        distances.push_back(min_distance);
    }
    int w = accumulate(distances.begin(), distances.end(), 0) / size;
    return w;
}