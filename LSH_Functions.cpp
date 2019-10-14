#include "lib/LSH_Functions.h"
#include "lib/Helper_Functions.h"

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
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    uniform_int_distribution<int> distribution (0, w);
    vector<int> Sj;

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < d; j++) {
            Sj.push_back(distribution(generator));
        }
        s->push_back(Sj);
        Sj.clear();
    }
}

void projections(vector<vector<int>>* a_projects, vector<vector<int>> x, vector<int>* s, int w, int d) {
    /* Ai = (Xi - Si) / W
     * Project every X to A in d-dimensional grid shifted by S, where every cell size = W */
    int ai;
    vector<int> a;
    for (int i = 0; i < x.size(); i++) {
        for (int dim = 0; dim < d; dim++) {
            ai = floor((double)(x[i][dim] - (*s)[dim]) / w) + w;
//            cout << " A Projects for " << i << " item is " << ai << endl;
            a.push_back(ai);
        }
        a_projects->push_back(a);
        a.clear();
    }
}

void compute_hash(vector<int>* H, vector<vector<int>> a, int d, int k, int w){
    /* we will compute K of hash functions for every point - item
     * vector H at the end will have a size of (dataset.size(), k) */
    /* TODO: we need to check for the size of every number -> has to be small, output G has to be 32bit */
    int m, M, h, term, term1, term2, fterm;
    m = w;
    M = pow(2, 32/k);
    for (int i = 0; i < a.size(); i++){
        h=0;
        for (int j = 0; j < d; j++){
            term = a[i][d-1-j]*pow(m,j);
            h += term % M;
//            term1 = a[i][d-1-j];
//            term2 = modpow<int>(m,j,M);
//            term = term1 * term2;
//            fterm = term;
//            h += fterm;
        }
//        cout << " Hash for " << i << " item is " << h << endl;
        H->push_back(h);
    }
}

void amplify_hash(vector<int>* amplified_g, vector<vector<int>>* hash_functions, int k){
    /* For every item it amplifies the hash from K dimensions to 1
     * g(x) = [h1(x)|h2(x)|h3(x)....|hk(x)] */
    int g;
    int concat_dist = 32/k;
    for (int i = 0; i < (*hash_functions)[0].size(); i++) {
        g=0;
        for (int j = 0; j < k; j++) {
            if (j == 0){
                g = ((*hash_functions)[j][i]);
            } else{
                g +=  g << concat_dist | ((*hash_functions)[j][i]);
            }
        }
        amplified_g->push_back(g);
        cout << "Amplified Hash for " << i << " item is " << (*amplified_g)[i] << endl;
    }
}

/* TODO: create brute force function to find the actual neighbors */