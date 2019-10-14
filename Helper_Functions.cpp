#include "lib/Helper_Functions.h"

using namespace std;

int Read_input_files(vector<vector<int>>* dataset, vector<vector<int>>* searchset, char* data_filename, char* query_filename) {
    string line;
    int number;
    vector<int> v;

    ifstream input_file(data_filename);
    if (!input_file) {
        cout << "Wrong input_file!" << endl;
        return -1;
    }
    ifstream query_file(query_filename);
    if (!query_file) {
        cout << "Wrong query_file!" << endl;
        return -1;
    }

    while (getline(input_file, line)) {
        stringstream ss(line);
        /* first number is the base */
        ss >> number;
        while (ss >> number) {
            v.push_back(number);
        }
        dataset->push_back(v);
        v.clear();
    }

    while (getline(query_file, line)) {
        stringstream ss(line);
        /* first number is the base */
        ss >> number;
        while (ss >> number) {
            v.push_back(number);
        }
        searchset->push_back(v);
        v.clear();
    }

    return 1;
}

int dist(vector<int>* P1, vector<int>* P2, int d, int Metric){
    /* Lk metric
     * for metric = 1 we have L1 metric
     * for metric = 2 we have L2 metric etc.
     * (default value = L1 Metric) -> Manhattan distance */
    int dist = 0;
    for (int dim = 0; dim < d; dim++)
        dist += pow(abs((*P1)[dim] - (*P2)[dim]),1);
    return pow(dist,1/Metric);
}

