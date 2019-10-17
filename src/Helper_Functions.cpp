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

int Read_curve_files(vector<vector<double*>>* dataset, char* data_filename) {
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

double dist(double* p, double* q, int Metric){
    /* Lk metric
     * for metric = 1 we have L1 metric
     * for metric = 2 we have L2 metric etc.
     * (default value = L1 Metric) -> Manhattan distance */
    double d1, d2;
    d1 = pow(p[0] - q[0], Metric);
    d2 = pow(p[1] - q[1], Metric);
    return pow(d1+d2,1/Metric);
}

double min(double x, double y, double z) {

}

void brute_force(vector<vector<int>>* dataset, vector<vector<int>>* searchset) {
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
        distances.push_back(min_distance);
        n_neighbors.push_back(n_neighbor);
        neighbors_file << "Item:" << i + 1 << ", Neighbor: " << n_neighbor << " | Distance: " << min_distance <<  endl;
    }
    neighbors_file.close();
}

void DTW(double *** c, vector<double*>* P, vector<double*>* Q) {
    /* Computing DTW
     * For m1 curves P and m2 curves Q, DTW computed in O(m1*m2) time and space by DP using recursion
     * Initialize c(1,1) = ||p1-q1||
     * if j > 1, then c(1,j) = c(1,j-1) + ||pi - qj||
     * if i > 1, then c(1,i) = c(i-1,1) + ||pi - qj||
     * if i > 1, j > 1, then c(i,j) = min{c(i-1,j), c(i-1,j-1), c(i,j-1)} + ||pi-qj|| */
    int m1, m2;
    int i,j;
    m1 = P->size();
    m2 = Q->size();
    (*c)[0][0] = dist((*P)[0], (*Q)[0], 2);
    if (j > 1 && i == 0) {
        (*c)[i][j] = (*c)[i][j - 1] + dist((*P)[i], (*Q)[j], 2);
    } else if (i > 1 && j == 0) {
        (*c)[i][j] = (*c)[i - 1][j] + dist((*P)[i], (*Q)[j], 2);
    } else {
        (*c)[i][j] = min((*c)[i-1][j], (*c)[i-1][j-1], (*c)[i][j-1]) + dist((*P)[i], (*Q)[j], 2);
    }
}