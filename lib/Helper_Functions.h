#include "Library.h"

using namespace std;

int Read_point_files(vector<vector<int>>*, vector<vector<int>>*, string, string );
int Read_curve_files(vector<vector<double*>>*, vector<vector<double*>>*, char*, char*);
template <typename Point> double dist(vector<Point>*, vector<Point>*, int, int=1);
double point_dist(double*, double*, int);
double DTW(vector<double*>*, vector<double*>*);
int modulo (int a, int b);
int moduloMultiplication(int, int, int);
long moduloPow(long, long, long);
void brute_force(vector<vector<int>>*, vector<vector<int>>*, vector<int>*, vector<double>*);
void curves_brute_force(vector<vector<double*>>*, vector<vector<double*>>*, vector<double>*, vector<double>*);
double* arg_min(double**, vector<double>*, double, int);
double min(double, double, double);

