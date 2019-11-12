#include "Initializers.h"
#include "Library.h"

using namespace std;

template <class Point>
void Random_Selection<Point>::init(vector<vector<Point>>* dataset) {

}

template <class Point>
string Random_Selection<Point>::get_name() {
    return this->name;
}

template <class Point>
void KMeans_plusplus<Point>::init(vector<vector<Point>>* dataset) {
    /* Choose a centroid uniformly at random;t←1 */
    for (int i = 0; i < this->K; i++) {

    }
    /* for all non-centroid point i=1,...,n−t, letD(i)←min distance to some centroid,
     * among t chosen centroids. */
    for (int i = 0; i < dataset->size(); i++) {

    }

    /* Choose new centroid: r chosen with probability proportional to D(r)^2 */
    for (int i = 0; i < this->K; i++) {

    }

    /* Go to (2) until t = k = given #centroids. */
}

template <class Point>
string KMeans_plusplus<Point>::get_name() {
    return this->name;
}

template class KMeans_plusplus<int>;
template class KMeans_plusplus<double>;
template class Random_Selection<int>;
template class Random_Selection<double>;