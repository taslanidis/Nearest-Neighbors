#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <bits/stdc++.h>               //for stringstream

#include "window.h"
#include "compute_s.h"
#include "compute_a.h"
#include "compute_H.h"
#include "compute_G.h"

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 3) {
      cout << "We need input_file AND query_file!" << endl;
      return -1;
    }
    int line_count = 0;
    int number_count = 0;
    string line;
    int number;

    ifstream input_file(argv[1]);
    if(!input_file){
        cout <<"Wrong input_file!" << endl;
        return -1;
    }
    ifstream query_file(argv[2]);
    if(!query_file){
        cout <<"Wrong query_file!" << endl;
        return -1;
    }
    vector<vector<int>> dataset;
    vector<vector<int>> searchset;
    while (getline(input_file, line)){
        vector<int> v;
        stringstream ss(line);
        ss >> number;                    //first number is the base
        while (ss >> number) {
            v.push_back(number);
        }
        dataset.push_back(v);
        v.clear();
        // if (line_count == 0){
        //    stringstream ss(line);
        //    while (ss >> number) {
        //      ++number_count;
        //    }
        // }
        // line_count++;
    }
    while (getline(query_file, line)){
        vector<int> v;
        stringstream ss(line);
        ss >> number;                    //first number is the base
        while (ss >> number) {
            v.push_back(number);
        }
        searchset.push_back(v);
        v.clear();
    }
    // cout << dataset[0].size() << endl;
    // cout << dataset.size() << endl;
    // cout << dataset.capacity() << endl;
    // cout << "Lines are: " << line_count << endl;
    // cout << "Words per line are: " << number_count << endl;

    int d = dataset[0].size();                           //d-dimensonal vectors

    int w = compute_window(dataset);

    vector<int> s;
    compute_s(&s, w, d);

    vector<vector<int>> a_projects;
    compute_a(&a_projects, dataset, &s, w, d);

    vector<int> H;
    int k = 4;
    compute_H(&H, a_projects, d, k, w);

    string gx = amplify(&H);

    cout << "G(x) = " << gx << endl;

    return 0;
}