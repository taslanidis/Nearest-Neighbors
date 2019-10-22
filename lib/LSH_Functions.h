#include "Library.h"

using namespace std;

template <typename Point> Point compute_window(vector<vector<Point>>*);
template <typename Point> void projections(vector<vector<int>>*, vector<vector<Point>>*, vector<double>*, int, int);

void generate_shifts(vector<vector<double>>*, int, int, int);
void compute_hash(vector<int>*, vector<vector<int>>*, int, int, int);
void amplify_hash(vector<int>*, vector<vector<int>>*, int);
