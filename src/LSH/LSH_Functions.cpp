#include "LSH_Functions.h"
#include "Helper_Functions.h"

using namespace std;

int compute_window(vector<vector<int>> dataset) {
    /* 1. take all points in dataset
     * 2. find their nearest neighbor using L1 metric
     * 3. average all distances between points
     * L1 = sum(|P1i - P2i|) for i in 0,d-1
     * Note: dist is a function scalable for all metrics(L1,L2,L3 etc.) */
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
                /* default is L1 metric, for Lk metric, add a 4th argument, k */
                L1 = dist(&P1, &P2, d);
                if (min_distance == -1)
                    min_distance = L1;
                if (L1 < min_distance)
                    min_distance = L1;
            }
            L1 = 0;
        }
        distances.push_back(min_distance);
    }
    int w = 10 * accumulate(distances.begin(), distances.end(), 0) / size;
    return w;
}

void generate_shifts(vector<vector<int>>* s, int w, int d, int k){
    /* Generate K * Si for every dimension
     * At the end, s will be a vector of size (k,d) */
    unsigned seed;

    uniform_int_distribution<int> distribution (0, w);
    vector<int> Sj;

    for (int i = 0; i < k; i++) {
        seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator (seed);
        for (int j = 1; j < d; j++) {
            Sj.push_back(distribution(generator));
        }
        s->push_back(Sj);
        Sj.clear();
    }
}

void projections(vector<vector<int>>* a_projects, vector<vector<int>>* x, vector<int>* s, int w, int d) {
    /* Ai = (Xi - Si) / W
     * Project every X to A in d-dimensional grid shifted by S, where every cell size = W */
    int ai;
    vector<int> a;
    for (int i = 0; i < x->size(); i++) {
        for (int dim = 1; dim < d; dim++) {
            ai = floor((double)((*x)[i][dim] - (*s)[dim - 1]) / w) + w; // TODO : used to be + w
            a.push_back(ai);
        }
        a_projects->push_back(a);
        a.clear();
    }
}

int power_of(int m, int j, int M) {
    int res = 1;
    m = m % M;
    while (j > 0)
    {
        if (j & 1)
            res = (res*m) % M;
        j = j>>1;
        m = (m*m) % M;
    }
    return res;
}

void compute_hash(vector<int>* H, vector<vector<int>> *a, int d, int k, int w){
    /* we will compute K of hash functions for every point - item
     * vector H at the end will have a size of (dataset.size(), k) */
    /* TODO: we need to check for the size of every number -> has to be small, output G has to be 32bit */
    int m, m1, m2, m3, m4, M, h, term;
    M = pow(2, 32/k);
    m1 = 2*16 % M;
    m2 = 2*8 % M;
    m3 = 2*4 % M;
    m4 = 2*2 % M;
    m = m1*m2*m3*m4*m4 % M;
    for (int i = 0; i < a->size(); i++){
        h=0;
        for (int j = 0; j < d - 1; j++){
            term = (*a)[i][d-1-j]*power_of(m,j,M); // TODO: warning check, previously had int j = 1
            h += term % M;  // moding with M to avoid overflow;
        }
        h = h % M;
        H->push_back(h);
    }
}

void amplify_hash(vector<int>* amplified_g, vector<vector<int>>* hash_functions, int k){
    /* For every item it amplifies the hash from K dimensions to 1
     * g(x) = [h1(x)|h2(x)|h3(x)....|hk(x)] */
    int g;
    int concat_dist = floor(32/k);
    /* for all points in dataset */
    for (int i = 0; i < (*hash_functions)[0].size(); i++) {
        g=0;
        for (int j = 0; j < k; j++) {
            if (j == 0) {
                g = (*hash_functions)[j][i];
            } else {
                //g +=  g << concat_dist | (*hash_functions)[j][i];
                g = g | ((*hash_functions)[j][i]); // different approach
            }
        }
        amplified_g->push_back(abs(g));
    }
}