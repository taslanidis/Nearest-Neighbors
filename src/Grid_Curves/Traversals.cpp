#include "Traversals.h"
#include "Helper_Functions.h"

void show_grid_lsh_usage(string name)
{
    cerr      << "Usage:   " << name << " -letter(s) <option(s)>\n"
              << "Options:\n"
              << "\t-d <input file>  Path to data file\n"
              << "\t-q <query file>  Path to search file\n"
              << "\t-k_vec <int>         Number of hi function for construction of g function\n"
              << "\t-L_grid <int>         Number of Grids\n"
              << "\t-o <output file> Path to file of results\n"
              << endl;
}

void show_grid_bhc_usage(string name)
{
    cerr      << "Usage:   " << name << " -letter(s) <option(s)>\n"
              << "Options:\n"
              << "\t-d <input file>  Path to data file\n"
              << "\t-q <query file>  Path to search file\n"
              << "\t-k_hypercube <int>         Dimension of Hypercube\n"
              << "\t-L_grid <int>         Number of Grids\n"
              << "\t-M <int>         Max allowed points to check\n"
              << "\t-probes <int>    Max allowed vertices to check\n"
              << "\t-o <output file> Path to file of results\n"
              << endl;
}

void shift_grid(vector<double>* orthogonal_grid, int delta, int d) {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    /* t uniformly in [0,d) */
    uniform_real_distribution<double> distribution (0, d);
    for (int i = 0; i < d; i++) {
        orthogonal_grid->push_back(distribution(generator));
    }
}

void hash_curve(vector<vector<double>>* hashed_curve, vector<double*>* curve, vector<double>* orthogonal_grid, double delta, int d) {
    double* pi;
    vector<double> pi_new;
    vector<double> pi_old;
    /* first index has the id and the dimensions of the curve */
    for (int i = 0; i < (*curve)[0][1]; i++) {
        pi = (*curve)[i];
        if (i == 0) {
            pi_new.push_back(pi[0]);
            pi_new.push_back(pi[1]);
            hashed_curve->push_back(pi_new);
        } else {
            pi_new = arg_min(&pi, orthogonal_grid, delta, d);
            /* remove consecutive duplicates pi' from the hashed_curve */
            if ((pi_new[0] != pi_old[0]) || (pi_new[1] != pi_old[1])) {
                hashed_curve->push_back(pi_new);
            }
        }
        pi_old = pi_new;
    }
}

void find_diagonal(vector<vector<int>>* v, int len1, int len2){
    double j = 0;
    double lamda = (double) len2 / len1;
    double stepi;
    if (len1>len2){
        stepi = (double) len2 / len1;
    } else{
        stepi = (double) len1 / len2;
    }
    for (double i = 0; i < len1; i+=stepi){
        j = lamda * i;
        vector<int> pair1;
        pair1.push_back(floor(i));
        pair1.push_back(floor(j));;
        vector<int> pair2;
        pair2.push_back(floor(i));
        pair2.push_back(floor(j));
        v->push_back(pair1);
        if((ceil(j) != floor(j)) && (ceil(j) < len2)) v->push_back(pair2);
    }
}

void Relevant_Traversals(vector<vector<vector<int>>>* Traversals, int len1, int len2) {
    vector<vector<int>> traversal;
    find_diagonal(&traversal, len1, len2);
    Traversals->push_back(traversal);
    /* find all +1 and -1 from the diagonal and push into vector */
}
