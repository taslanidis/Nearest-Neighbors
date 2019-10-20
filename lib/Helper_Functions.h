#include "Library.h"

using namespace std;

int Read_point_files(vector<vector<int>>*, vector<vector<int>>*, char*, char* );
int Read_curve_files(vector<vector<double*>>*, vector<vector<double*>>*, char*, char*);
int dist(vector<int>*, vector<int>*, int, int=1);
double point_dist(double*, double*, int);
int modulo (int a, int b);
long long moduloMultiplication(long long a, long long b, long long mod);
void brute_force(vector<vector<int>>*, vector<vector<int>>*, vector<int>*, vector<double>*);
double* arg_min(double**, vector<int>*, double, int);