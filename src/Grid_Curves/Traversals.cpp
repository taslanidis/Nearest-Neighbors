#include "Traversals.h"
#include "Helper_Functions.h"

void shift_grid(vector<double>* orthogonal_grid, int delta, int d) {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    /* t uniformly in [0,d) */
    uniform_real_distribution<double> distribution (0, d);
    for (int i = 0; i < d; i++) {
        orthogonal_grid->push_back(distribution(generator));
    }
}

void hash_curve(vector<double*>* hashed_curve, vector<double*>* curve, vector<double>* orthogonal_grid, double delta, int d) {
    double* pi;
    double* pi_new;
    double* pi_old;
    /* first index has the id and the dimensions of the curve */
    for (int i = 0; i < (*curve)[0][1]; i++) {
        pi = (*curve)[i];
        if (i == 0) {
            pi_new = pi;
            pi_old = pi_new;
        } else {
            pi_new = arg_min(&pi, orthogonal_grid, delta, d);
        }
        /* remove consecutive duplicates pi' from the hashed_curve
        if (pi_new[0] != pi_old[0] && pi_new[1] != pi_old[1]) */
        hashed_curve->push_back(pi_new);
        pi_old = pi_new;
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