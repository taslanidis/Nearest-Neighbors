#include "Library.h"
#include "Cluster.h"

using namespace std;

template <class Point>
Cluster <Point>::Cluster(int K, string Initializer, string Assigner) {
    /* k for k means */
    this->K = K;
    /* Initializer */
    if (Initializer == "Random Selection") {
        this->initializer = new Random_Selection<Point>(K);
        cout << initializer->get_name() << endl;
    } else if (Initializer == "K-Means++") {
        this->initializer = new KMeans_plusplus<Point>(K);
        cout << initializer->get_name() << endl;
    } else {
        cerr << "Unknown Initializer";
    }
    /* Assigner */
    if (Assigner == "Lloyd's Assignment") {
        this->assigner = new Lloyd_assignment<Point>();
        cout << assigner->get_name() << endl;
    } else if (Assigner == "Inverse Assignment") {
        this->assigner = new Inverse_assignment<Point>();
        cout << assigner->get_name() << endl;
    } else {
        cerr << "Unknown Assigner";
    }
}

template <class Point>
void Cluster <Point>::fit(vector<vector<Point>>* dataset) {
    /* initialization */
    cout << '\t' << "Initializer call ..." << endl;
    this->centroids = initializer->init(dataset);
    /* assignment */
    cout << '\t' << "Assigner call ..." << endl;
    assigner->assign(dataset, this->centroids);
    /* update */
}

template <class Point>
Cluster <Point>::~Cluster(){
    delete (this->initializer);
}

template class Cluster<int>;
template class Cluster<double>;