#include "Traversals.h"
#include "Helper_Functions.h"

double min(double x, double y, double z) {
    /* get the min out of 3 real numbers */
    double temp = (x < y) ? x : y;
    return (z < temp) ? z : temp;
}

void Traversal(double*** c, vector<double*>* P, vector<double*>* Q) {
    /* Initialize c(1,1) = ||p1-q1||
     * * if j > 1, then c(1,j) = c(1,j-1) + ||pi - qj||
    * if i > 1, then c(1,i) = c(i-1,1) + ||pi - qj||
    * if i > 1, j > 1, then c(i,j) = min{c(i-1,j), c(i-1,j-1), c(i,j-1)} + ||pi-qj|| */
    /* TODO: make it with dynamic programming, not exhaustive */
    int m1 = P->size();
    int m2 = Q->size();
    (*c)[0][0] = point_dist((*P)[0], (*Q)[0], 2);
    for (int i = 1; i < m1; i++) {
        for (int j = 1; j < m2; j++) {
            if (j > 0 && i == 0) {
                (*c)[i][j] = (*c)[i][j - 1] + point_dist((*P)[i], (*Q)[j], 2);
            } else if (i > 0 && j == 0) {
                (*c)[i][j] = (*c)[i - 1][j] + point_dist((*P)[i], (*Q)[j], 2);
            } else {
                (*c)[i][j] = min((*c)[i - 1][j], (*c)[i - 1][j - 1], (*c)[i][j - 1]) + point_dist((*P)[i], (*Q)[j], 2);
            }
        }
    }
}

void DTW(double *** c, vector<double*>* P, vector<double*>* Q) {
    /* Computing DTW
     * For m1 curves P and m2 curves Q, DTW computed in O(m1*m2) time and space by DP using recursion */
    Traversal(c, P, Q);
    // find min traversal from c
}