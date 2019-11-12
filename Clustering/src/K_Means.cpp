#include "Library.h"
#include "K_Means.h"

using namespace std;

template <class Point>
K_Means <Point>::K_Means(int K, string Initializer) {
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
K_Means <Point>::~K_Means(){
    delete (this->initializer);
}

template class K_Means<int>;
template class K_Means<double>;