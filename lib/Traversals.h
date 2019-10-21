#include "Library.h"

using namespace std;

double DTW(double ***, vector<double*>*, vector<double*>*);
double min(double, double, double);
void shift_grid(vector<int>*, int, int);
void hash_curve(vector<double*>*, vector<double*>*, vector<int>*, double, int);
void Relevant_Traversals(vector<int*>*, int, int, int, int);