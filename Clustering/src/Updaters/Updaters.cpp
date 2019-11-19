#include "Updaters.h"
#include "Library.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
void PAM<Point>::update(vector<vector<Point>>* dataset, vector<int>** clusters, vector<int>* centroids) {
    /* minimize Sum(dist(i,t)) over all objects t in cluster C */
    /* todo: optimizations for speed */
    int t;
    double sum, min;
    for (int i = 0; i < this->get_K(); i++) {
        cout << "PAM for cluster <" << i << "> " << endl;
        min = DBL_MAX;
        for (int j = 0; j < clusters[i]->size(); j++) {
            /* find medoid t to minimize the distances in this cluster */
            sum = 0.0;
            for (int l = 0; l < clusters[i]->size(); l++) {
                if (j == l) continue;
                sum += dist(&(*dataset)[(*clusters[i])[j]], &(*dataset)[(*clusters[i])[l]]);
            }
            /* find min and the id of min, make it centroid for this cluster */
            if (sum < min) {
                min = sum;
                t = (*clusters[i])[j];
            }
        }
        /* new centroid for this cluster */
        centroids->at(i) = t;
    }
    return;
}

template <class Point>
string PAM<Point>::get_name() {
    return this->name;
}

template <class Point>
void MV_DTW<Point>::update(vector<vector<Point>>* dataset, vector<int>** clusters, vector<int>* centroids) {
    return;
}

template <class Point>
string MV_DTW<Point>::get_name() {
    return this->name;
}

template class PAM<int>;
template class PAM<double*>;
template class MV_DTW<int>;
template class MV_DTW<double*>;