#include "Helper_Functions.h"

using namespace std;

void show_cluster_usage(string name)
{
    cerr      << "Usage:   " << name << " -letter(s) <option(s)>\n"
              << "Options:\n"
              << "\t-i <input file>             Path to input.dat\n"
              << "\t-c <configuration file>     Path to cluster.conf\n"
              << "\t-o <output file>            Path to file of results\n"
              << "\t-complete <optional>        Thorough Results\n"
              << endl;
}

int Read_input_file(string input){
    string line;
    ifstream input_file(data_filename);
    getline(input_file, line);
    if (line == "vectors"){
        return 1;
    }else if (line == "curves"){
        return 2;
    }else{
        return 0;
    }
}

int Read_point_files(vector<vector<int>>* dataset, vector<vector<int>>* searchset, string data_filename, string query_filename) {
    /* Create dataset and searchset vectors from files */

    string line;
    int id;
    int number;
    vector<int> v;

    ifstream input_file(data_filename);
    ifstream query_file(query_filename);

    id = 0;
    while (getline(input_file, line)) {
        stringstream ss(line);
        /* discard id */
        ss >> number;
        v.push_back(id);
        while (ss >> number) {
            v.push_back(number);
        }
        dataset->push_back(v);
        v.clear();
        id++;
    }

    id = -1;
    while (getline(query_file, line)) {
        stringstream ss(line);
        if(id == -1){
            /* Radius: */
            ss >> line;
            ss >> (*R);
            id++;
        }else {
            /* discard id */
            ss >> number;
            v.push_back(id);
            while (ss >> number) {
                v.push_back(number);
            }
            searchset->push_back(v);
            v.clear();
            id++;
        }
    }

    if (dataset->size() == 0) {
        cerr << "Error: dataset file is empty!" << endl;
        return -1;
    }
    /* <= 1 because first line is radius */
    if (searchset->size() <= 1) {
        cerr << "Error: searchset file is empty!" << endl;
        return -1;
    }

    return 1;
}

int Read_curve_files(vector<vector<double*>>* dataset, vector<vector<double*>>* searchset, string data_filename, string query_filename) {
    string line;
    int id, trash;
    char bracket, comma;
    double number;
    double* point;
    vector<double*> v;

    ifstream input_file(data_filename);
    ifstream query_file(query_filename);

    id = 0;
    while (getline(input_file, line)) {
        stringstream ss(line);
        /* id */
        ss >> trash;
        point = new double [2];
        point[0] = id;
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
        id++;
    }

    id = 0;
    while (getline(query_file, line)) {
        stringstream ss(line);
        /* id */
        ss >> trash;
        point = new double [2];
        point[0] = id;
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
        id++;
    }

    if (dataset->size() == 0) {
        cerr << "Error: dataset file is empty!" << endl;
        return -1;
    }
    if (searchset->size() == 0) {
        cerr << "Error: searchset file is empty!" << endl;
        return -1;
    }
    return 1;
}