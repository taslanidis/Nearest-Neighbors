#include "Helper_Functions.h"
#include "HashTable.h"

using namespace std;

int Read_point_files(vector<vector<int>>* dataset, vector<vector<int>>* searchset, char* data_filename, char* query_filename) {
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
        while (ss >> number) {
            v.push_back(number);
        }
        dataset->push_back(v);
        v.clear();
    }

    while (getline(query_file, line)) {
        stringstream ss(line);
        while (ss >> number) {
            v.push_back(number);
        }
        searchset->push_back(v);
        v.clear();
    }

    return 1;
}

int Read_curve_files(vector<vector<double*>>* dataset, vector<vector<double*>>* searchset, char* data_filename, char* query_filename) {
    string line;
    char bracket, comma;
    double number;
    double* point;
    vector<double*> v;

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
        /* id */
        point = new double [2];
        ss >> point[0];
        /* length */
        ss >> point[1];
        v.push_back(point);
        /* for all data points */
        while (ss >> bracket) {
            point = new double [2];
            ss >> point[0];
            ss >> comma;
            ss >> point[1];
            ss >> bracket;
            v.push_back(point);
        }
        dataset->push_back(v);
        v.clear();
    }

    while (getline(query_file, line)) {
        stringstream ss(line);
        /* id */
        point = new double [2];
        ss >> point[0];
        /* length */
        ss >> point[1];
        v.push_back(point);
        /* for all data points */
        while (ss >> bracket) {
            point = new double [2];
            ss >> point[0];
            ss >> comma;
            ss >> point[1];
            ss >> bracket;
            v.push_back(point);
        }
        searchset->push_back(v);
        v.clear();
    }

    return 1;
}

int dist(vector<int>* P1, vector<int>* P2, int d, int Metric) {
    /* Lk metric
     * for metric = 1 we have L1 metric
     * for metric = 2 we have L2 metric etc.
     * (default value = L1 Metric) -> Manhattan distance */
    int dist = 0;
    for (int dim = 1; dim < d; dim++)
        dist += pow(abs((*P1)[dim] - (*P2)[dim]),Metric);
    return pow(dist,1/Metric);
}

double point_dist(double* p, double* q, int Metric){
    /* Lk metric
     * for metric = 1 we have L1 metric
     * for metric = 2 we have L2 metric etc.
     * (default value = L1 Metric) -> Manhattan distance */
    double d1, d2;
    d1 = pow(p[0] - q[0], Metric);
    d2 = pow(p[1] - q[1], Metric);
    return pow(d1+d2,1/Metric);
}

void brute_force(vector<vector<int>>* dataset, vector<vector<int>>* searchset, vector<int>* TrueDistances, vector<double>* TrueTimes) {
    /* vectors init */
    vector<int> n_neighbors;
    vector<int> distances;
    vector<int> P1;
    vector<int> P2;
    /* variable declaration */
    int L1, min_distance, n_neighbor;
    int d_size = dataset->size();
    int s_size = searchset->size();
    int d = (*dataset)[0].size();
    /* open file to write results */
    ofstream neighbors_file;
    neighbors_file.open ("./output/nneighbors_bf.txt");
    for (int i = 0; i < s_size; i++) {
        P1 = (*searchset)[i];
        L1 = 0;
        min_distance = -1;
        auto start = chrono::high_resolution_clock::now();
        for (int j = 1; j < d_size; j++) {
            P2 = (*dataset)[j];
            /* default is L1 metric, for Lk metric, add a 4th argument, k */
            L1 = dist(&P1, &P2, d);
            if (min_distance == -1) {
                min_distance = L1;
                n_neighbor = j;
            }
            if (L1 < min_distance) {
                min_distance = L1;
                n_neighbor = j;
            }
        }
        auto finish = chrono::high_resolution_clock::now();
        auto elapsed = finish - start;
        double time_elapsed = chrono::duration<double>(elapsed).count();
        TrueDistances->push_back(min_distance);
        TrueTimes->push_back(time_elapsed);
        distances.push_back(min_distance);
        n_neighbors.push_back(n_neighbor);
        neighbors_file << "Item:" << setw(floor(log10(s_size) + 1)) << setfill('0') << i + 1 << ", Neighbor: " << setw(floor(log10(d_size) + 1)) << setfill('0') << n_neighbor + 1 << " | Distance: " << setw(4) << setfill('0') << min_distance <<  " | Duration: " << time_elapsed << endl;
    }
    neighbors_file.close();
}

double* arg_min(double** pi, vector<int>* orthogonal_grid, double delta, int d) {
    double min, norm;
    int q;
    double* argmin = new double[d];
    /* Point is to minimize the ||pi-q|| for all q
     * Solve this for zero ->
     * delta * ai + (*orthogonal_grid)[i]) - (*pi)[0] = 0*/
    for (int i = 0; i < d; i++) {
        argmin[i] = abs((*orthogonal_grid)[i] - (*pi)[i]) / delta + (*orthogonal_grid)[i];
    }
    return argmin;
}