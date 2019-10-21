#include "Traversals.h"
#include "Helper_Functions.h"

double min(double x, double y, double z) {
    /* get the min out of 3 real numbers */
    double temp = (x < y) ? x : y;
    return (z < temp) ? z : temp;
}

double DTW(double *** c, vector<double*>* P, vector<double*>* Q) {
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
    return (*c)[m1-1][m2-1];
}

void shift_grid(vector<int>* orthogonal_grid, int delta, int d) {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    /* t uniformly in [0,d) */
    uniform_int_distribution<int> distribution (0, d);
    for (int i = 0; i < d; i++) {
        orthogonal_grid->push_back(distribution(generator));
    }
}

void hash_curve(vector<double*>* hashed_curve, vector<double*>* curve, vector<int>* orthogonal_grid, double delta, int d) {
    double* pi;
    double* pi_new;
    /* first index has the id and the dimensions of the curve */
    for (int i = 0; i < (*curve)[0][1]; i++) {
        pi = (*curve)[i];
        if (i == 0) {
            pi_new = pi;
        } else {
            pi_new = arg_min(&pi, orthogonal_grid, delta, d);
        }
        /* TODO: remove consecutive duplicates pi' from the hashed_curve */
        hashed_curve->push_back(pi_new);
    }
}

/* Find the ORANGE boxes in the array
 * we keep the indexes of the points of every curve
 * TODO: when dimensions are very different this doesnt go well */
void Relevant_Traversals(vector<int*>* Traversals, int len1, int len2, int i, int j) {
    int * pair = new int [2];
    pair[0] = i;
    pair[1] = j;
    Traversals->push_back(pair);
    if (i == len1 && j == len2 ) {
        return;
    } else if ( j == len2) {
        Relevant_Traversals(Traversals, len1, len2, i + 1, j);
    } else if ( i == len1 ) {
        Relevant_Traversals(Traversals, len1, len2, i, j + 1);
    } else if (i == j) {
        Relevant_Traversals(Traversals, len1, len2, i + 1, j);
        Relevant_Traversals(Traversals, len1, len2, i, j + 1);
    } else if (i < j) {
        Relevant_Traversals(Traversals, len1, len2, i+1, j);
    } else {
        Relevant_Traversals(Traversals, len1, len2, i, j + 1);
    }
}