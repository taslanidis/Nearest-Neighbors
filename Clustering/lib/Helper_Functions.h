#include "Library.h"

using namespace std;

template <typename Point> double min_distance(int, vector<int>*, vector<vector<Point>>*);

void show_cluster_usage(string name);

int Read_input_file(string input);
int Read_files(vector<vector<int>>*, int*, string, string );
int Read_files(vector<vector<double*>>*, int*, string, string);
void normalize(vector<double>*);
double Sum(int, int, vector<double>*, int);
