#include "Updaters.h"
#include "Library.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
void PAM<Point>::update() {
    return;
}

template <class Point>
string PAM<Point>::get_name() {
    return this->name;
}

template <class Point>
void MV_DTW<Point>::update() {

}

template <class Point>
string MV_DTW<Point>::get_name() {
    return this->name;
}

template class PAM<int>;
template class PAM<double*>;
template class MV_DTW<int>;
template class MV_DTW<double*>;