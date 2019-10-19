#include "BHC_Functions.h"

using namespace std;

void fill_dictionary (map<int,int>* dictionary,vector<vector<int>> amplified_g) {
    default_random_engine generator;
    uniform_int_distribution<int> distribution(0,1);
    map<int,int>::iterator map_it;
//    cout << amplified_g.size() << endl;
//    cout << amplified_g[0].size() << endl;
    for (int i = 0; i < amplified_g.size(); i++) {
        for (int j=0; j < amplified_g[0].size(); j++) {
            map_it = dictionary->find(amplified_g[i][j]);
            if (map_it == dictionary->end()) {
                dictionary->insert(pair<int, int>(amplified_g[i][j], distribution(generator)));
            }
        }
    }
}

int calculate_vertex (vector<vector<int>> amplified_g, map<int,int> dictionary, int index) {
    int vertex = 0;
    map<int,int>::iterator map_it;
    for (int i=0; i < amplified_g.size(); i++){
        map_it = dictionary.find(amplified_g[i][index]);
        if (map_it != dictionary.end()) {
//            vertex += to_string(map_it->second);
            vertex = vertex << 1 | map_it->second;
        }
    }
    return vertex;
}

bool check_compatibility(vector<vector<int>>* query_amplified_g, vector<vector<int>>* data_amplified_g, int query_index, int data_index) {
    bool check = true;
    for(int i = 0; i < (*query_amplified_g).size(); i++){
        if((*query_amplified_g)[i][query_index] != (*data_amplified_g)[i][data_index]){
//            cout << "Query = " << (*query_amplified_g)[i][query_index] << endl;
//            cout << "Data = " << (*data_amplified_g)[i][data_index] << endl;
//            getchar();
            check = false;
            break;
        }
    }
    return check;
}
