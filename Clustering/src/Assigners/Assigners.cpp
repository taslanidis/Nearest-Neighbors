#include "Assigners.h"
#include "Library.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
vector<vector<int>>* Lloyd_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<int>* centroids) {
    int num_of_centroids = centroids->size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size();
    vector<vector<int>>* clusters = new vector<vector<int>>;
    int centroid;
    double min_dist, curr_dist;
    if(typeid(Point) == typeid(int)) {
        for (int i = 0; i < data_size; i++) {
            if (find(centroids->begin(), centroids->end(), i) != centroids->end()) continue;
            centroid = -1;
            min_dist = DBL_MAX;
            for (int j = 0; j < num_of_centroids; j++) {
//                curr_dist = dist(&dataset[(*centroids)[j]], &dataset[i], dimension);              //todo: causes error to makefile
                if (curr_dist < min_dist) {
                    min_dist = curr_dist;
                    centroid = j;
                }
            }
            (*clusters)[centroid].push_back(i);
        }
    }else if(typeid(Point) == typeid(double)){

    }else{
        cerr << "Wrong Point type in Lloyd's Assignment" << endl;
    }
    return clusters;
}

template <class Point>
string Lloyd_assignment<Point>::get_name() {
    return this->name;
}

template <class Point>
vector<vector<int>>* Inverse_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<int>* centroids) {

}

template <class Point>
string Inverse_assignment<Point>::get_name() {
    return this->name;
}

template class Lloyd_assignment<int>;
template class Lloyd_assignment<double*>;
template class Inverse_assignment<int>;
template class Inverse_assignment<double*>;