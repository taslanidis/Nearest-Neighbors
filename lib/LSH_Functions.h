#include "Library.h"

using namespace std;

int compute_window(vector<vector<int>>);
void projections(vector<vector<int>>*, vector<vector<int>>*, vector<int>*, int, int);
void generate_shifts(vector<vector<int>>*, int, int, int);
void compute_hash(vector<int>*, vector<vector<int>>*, int, int, int);
void amplify_hash(vector<int>*, vector<vector<int>>*, int);