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

    int* distances = (int*)malloc(dataset.size() * sizeof(int));
    for (int i = 0; i < dataset.size(); i++) {
        vector<int> P1 = dataset[i];
        vector<int> P2;
        int L1 = 0;
        int min_distance = -1;
        for (int j = 0; j < dataset.size(); j++) {
            if (i != j) {
                P2 = dataset[j];
                for (int dim = 0; dim < dataset[i].size(); dim++)
                    L1 += abs(P1[dim] - P2[dim]);
                if (min_distance == -1)
                    min_distance = L1;
                if (L1 < min_distance)
                    min_distance = L1;
            }
        }
        distances[i] = min_distance;
    }
    int w = accumulate(distances, distances + dataset.size(), 0)/dataset.size();
    free(distances);

    return w;
}