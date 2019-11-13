#include "Library.h"
#include "Cluster.h"

using namespace std;

template <class Point>
Cluster <Point>::Cluster(int K, string Initializer) {
    /* k for k means */
    this->K = K;
    /* initializer */
    if (Initializer == "Random Selection") {
        this->initializer = new Random_Selection<Point>(K);
        cout << initializer->get_name() << endl;
    } else if (Initializer == "K-Means++") {
        this->initializer = new KMeans_plusplus<Point>(K);
        cout << initializer->get_name() << endl;
    } else {
        cerr << "Unknown initializer";
    }
}

template <class Point>
void Cluster <Point>::fit(vector<vector<Point>>* dataset) {
    /* initialization */
    cout << '\t' << "Initializer call ..." << endl;
    initializer->init(dataset);
    /* assignment */
    cout << '\t' << "Assignment ..." << endl;
    /* update */
}

template <class Point>
Cluster <Point>::~Cluster(){
    delete (this->initializer);
}

template class Cluster<int>;
template class Cluster<double>;