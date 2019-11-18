#include "Initializers.h"
#include "Library.h"
#include "Helper_Functions.h"
#include <algorithm>

using namespace std;

template <class Point>
vector<int>* Random_Selection<Point>::init(vector<vector<Point>>* dataset) {
    vector<int>* centroids = new vector<int>;
    cout << '\t' << "Initializing with Random Selection" << endl;
    return centroids;
}

template <class Point>
string Random_Selection<Point>::get_name() {
    return this->name;
}

template <class Point>
vector<int>* KMeans_plusplus<Point>::init(vector<vector<Point>>* dataset) {
    cout << '\t' << "Initializing with K-Means++" << endl;
    /* # of centroids */
    int t = 1;
    /* # of data */
    int n = dataset->size();
    /* ids of points in dataset */
    vector<int>* centroids = new vector<int>;
    /* array with min distances of points from centroids */
    vector<double> D;
    /* vector for partial sums */
    vector<double> P;
    /* uniform distribution */
    unsigned seed;
    uniform_int_distribution<int> distribution (0, dataset->size()-1);
    seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    /* Step 1: Choose a centroid uniformly at random; t←1 */
    centroids->push_back(distribution(generator));
    do {
        /* Step 2: for all non-centroid point i=1,...,n−t, letD(i)←min distance to some centroid,
         * among t chosen centroids. */
        for (int i = 0; i < n; i++) {
            // for all non centroids
            if (find(centroids->begin(), centroids->end(), i) != centroids->end()) {
                D.push_back(min_distance(i, centroids, dataset));
            }else {
                D.push_back(0);
            }
        }

        /* Normalization */
        normalize(&D);

        /* Step 3: Choose new centroid: r chosen with probability proportional to D(r)^2 */
        for (int r = 1; r <= n - t; r++) {
            P.push_back(Sum(1, r, &D, 2));
        }

        /* sort vector P */
        sort(P.begin(), P.end());

        /* uniform distribution */
        unsigned seed;
        uniform_real_distribution<double> distribution(0, P[n - t - 1]);
        seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        /* Step 1: Choose a centroid uniformly at random; t←1 */
        double x = distribution(generator);

        int r = 0;
        for (r = 0; r < P.size(); r++) {
            if (x >= P[r]) break;
        }
        t++;
        /* TODO: careful, r is the index of the sorted data, i need to keep index of the initial ones */
        centroids->push_back(r);

        /* Step 4: Go to (2) until t = k = given #centroids. */
        vector<double>().swap(D);
        vector<double>().swap(P);
    }
    while (t < this->get_K());
    return centroids;
}

template <class Point>
string KMeans_plusplus<Point>::get_name() {
    return this->name;
}

template class KMeans_plusplus<int>;
template class KMeans_plusplus<double*>;
template class Random_Selection<int>;
template class Random_Selection<double*>;