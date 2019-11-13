#include "Assigners.h"
#include "Library.h"

using namespace std;

template <class Point>
void Lloyd_assignment<Point>::assign(vector<vector<Point>>* dataset) {

}

template <class Point>
string Lloyd_assignment<Point>::get_name() {
    return this->name;
}

template <class Point>
void Inverse_assignment<Point>::assign(vector<vector<Point>>* dataset) {

}

template <class Point>
string Inverse_assignment<Point>::get_name() {
    return this->name;
}

template class Lloyd_assignment<int>;
template class Lloyd_assignment<double>;
template class Inverse_assignment<int>;
template class Inverse_assignment<double>;