#include "Initializers.h"
#include "Library.h"

using namespace std;

void Random_Selection::init() {

}

string Random_Selection::get_name() {
    return this->name;
}

void KMeans_plusplus::init() {
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

string KMeans_plusplus::get_name() {
    return this->name;
}
