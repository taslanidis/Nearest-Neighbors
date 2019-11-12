#include "Library.h"

using namespace std;

void show_cluster_usage(string name);

int Read_input_file(string input);
int Read_point_files(vector<vector<int>>*, vector<vector<int>>*, string, string );
int Read_curve_files(vector<vector<double*>>*, vector<vector<double*>>*, string, string);