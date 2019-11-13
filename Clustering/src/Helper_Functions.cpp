#include "Helper_Functions.h"

using namespace std;

template double min_distance<int>(int, vector<int>*, vector<vector<int>>*);
template double min_distance<double>(int, vector<int>*, vector<vector<double>>*);

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
    ifstream input_file(input);
    getline(input_file, line);
    if (!line.empty() && line[line.size() - 1] == '\r')
        line.erase(line.size() - 1);
    if (!line.empty() && line[line.size() - 1] == '\n')
        line.erase(line.size() - 1);
    if (line == "vectors"){
        return 1;
    }else if (line == "curves"){
        return 2;
    }else{
        return 0;
    }
}

int Read_files(vector<vector<int>>* cluster_data, int* cluster_config, string input_file_name, string config_file_name) {

    string line;
    int id;
    int number;
    vector<int> v;

    ifstream input_file(input_file_name);
    ifstream config_file(config_file_name);

    id = 0;
    /* discard "vectors" */
    getline(input_file, line);
    while (getline(input_file, line)) {
        stringstream ss(line);
        /* discard id */
        ss >> number;
        v.push_back(id);
        while (ss >> number) {
            v.push_back(number);
        }
        cluster_data->push_back(v);
        v.clear();
        id++;
    }

    id = 0;
    while (getline(config_file, line)) {
        stringstream ss(line);
        while (line != ":"){
            ss >> line;
        }
        ss >> number;
        cluster_config[id] = number;
        id++;
    }

    if (cluster_data->size() == 0) {
        cerr << "Error: input.dat file is empty!" << endl;
        return -1;
    }

    return 1;
}

int Read_files(vector<vector<double*>>* cluster_data, int* cluster_config, string input_file_name, string config_file_name) {

    string line;
    int id, trash;
    char bracket, comma;
    double number;
    double* point;
    vector<double*> v;

    ifstream input_file(input_file_name);
    ifstream config_file(config_file_name);

    id = 0;
    /* discard "curves" */
    getline(input_file, line);
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
        cluster_data->push_back(v);
        v.clear();
        id++;
    }

    id = 0;
    while (getline(config_file, line)) {
        stringstream ss(line);
        while (line != ":"){
            ss >> line;
        }
        ss >> number;
        cluster_config[id] = number;
        id++;
    }

    if (cluster_data->size() == 0) {
        cerr << "Error: input.dat file is empty!" << endl;
        return -1;
    }
    return 1;
}

void normalize(vector<double>* D) {
    return;
}

double Sum(int start, int end, vector<double>* D, int power) {
    return 0.0;
}

template <typename Point>
double min_distance(int index, vector<int>* centroids, vector<vector<Point>>* dataset) {
    return 0.0;
}

