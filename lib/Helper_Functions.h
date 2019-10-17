#include "Library.h"

using namespace std;

int Read_point_files(vector<vector<int>>*, vector<vector<int>>*, char*, char* );
int Read_curve_files(vector<vector<double*>>*, char*);
int dist(vector<int>*, vector<int>*, int, int=1);
double dist(double*, double*, int=1);
void brute_force(vector<vector<int>>*, vector<vector<int>>*);
void DTW(double ***, vector<vector<double*>>*, vector<vector<double*>>*);