#include "lib/Library.h"
#include "lib/LSH_Functions.h"
#include "lib/BHC_Functions.h"
#include "lib/Helper_Functions.h"
#include "lib/HashTable.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "We need input_file AND query_file!" << endl;
        return -1;
    }
    /* variable declaration | k = 4 default value */
    int error_code, k = 4, dim = 3, M = 10, probes = 2;                       // k is the number of hi concatenated to form g - dim is number of hypercube's vertices
    /* vectors for the data and query points */
    vector<vector<int>> dataset;
    vector<vector<int>> searchset;

    /* read data set and query set and load them in vectors */
    error_code = Read_input_files(&dataset, &searchset, argv[1], argv[2]);
    if (error_code == -1) return -1;

    /* do brute force to find actual NNs */
#ifdef BRUTE_FORCE
    brute_force(&dataset, &searchset);
#endif

    /* d-dimensional vectors */
    int d = dataset[0].size();
    /* compute window for all hash tables (try *4 or *10) */
    //int w = 4*compute_window(dataset);
    int w = 36;
    /* vector containing (k,d) shifts */
    vector<vector<int>> s;
    /* H of size (k, dataset.size()) */
    vector<vector<int>> hash_functions;
    /* projections of data */
    vector<vector<int>> a_projects;
    /* internal vector H for pushing */
    vector<int> H;
    /* amplified hash for dataset*/
    vector<vector<int>> data_amplified_g;
    /* amplified hash for searchset */
    vector<vector<int>> query_amplified_g;
    /* map to assign every g function to 0 or 1 */
    map <int, int> dictionary;
    /* loop for L, to create L amplified functions g */
    for (int l = 0; l < dim; l++) {
        /* generate the random shifts */
        generate_shifts(&s, w, d, k);

        /* ----------------------- DATA SET -------------------------------*/
        /* loop for K */
        for (int i = 0; i < k; i++) {
            a_projects.clear();
            projections(&a_projects, dataset, &(s[i]), w, d);
            H.clear();
            compute_hash(&H, a_projects, d, k, w);
            hash_functions.push_back(H);
        } /* end for */
        /* compute the amplified hashes for every data item */
        amplify_hash(&data_amplified_g, &hash_functions, k);

        /* ------------------------ SEARCH SET ----------------------------*/
        /* loop for K */
        for (int i = 0; i < k; i++) {
            a_projects.clear();
            projections(&a_projects, searchset, &(s[i]), w, d);
            H.clear();
            compute_hash(&H, a_projects, d, k, w);
            hash_functions.push_back(H);
        } /* end for */
        /* compute the amplified hashes for every query item */
        amplify_hash(&query_amplified_g, &hash_functions, k);
    }

    fill_dictionary(&dictionary, data_amplified_g);
//    for(auto& x: dictionary){
//        cout << x.first << ", " << x.second << endl;
//    }

    int hypercube_size = pow(2,dim);
    int vertex = 0;
    HashTable *MyHashTable;
    MyHashTable = new HashTable(hypercube_size);
    for (int i = 0; i < dataset.size(); i++) {
        vertex = calculate_vertex(data_amplified_g, dictionary, i);
//        cout << i << " : " << vertex << endl;
        MyHashTable->Insert_to_Vertex(vertex, dataset[i]);
    }

    vector<vector<vector<int>>> ANNi;
    for (int i = 0; i < searchset.size(); i++) {
        vertex = calculate_vertex(query_amplified_g, dictionary, i);
//        cout << i << " : " << vertex << endl;
        ANNi.push_back(*MyHashTable->Search_Neighbors(vertex));
    }

    int distance = 0;
    int *min_distance = new int [searchset.size()];
    int *nearest_neighbor = new int [searchset.size()];

    /* initialize arrays */
    for (int i = 0; i < searchset.size(); i++) {
        min_distance[i] = INT_MAX;
        nearest_neighbor[i] = -1;
    }

    for (int q = 0; q < searchset.size(); q++) {
        /* for every vector in the same bucket (max M calculations) */
        for (int j = 0; j < ANNi[q].size() && j < M; j++) {
            /* TODO: I have to check for same g(x) also */
            distance = dist(&ANNi[q][j], &searchset[q], d);
            if (distance < min_distance[q]) {
                min_distance[q] = distance;
                nearest_neighbor[q] = ANNi[q][j][0];
            }
        }
    }

    /* print results */
    /* open file to write results */
    ofstream neighbors_file;
    neighbors_file.open ("nneighbors_bhc.txt");
    for (int i = 0; i < searchset.size(); i++) {
        neighbors_file << "Item: " << i + 1 << ", Neighbor: " << nearest_neighbor[i] << " | Distance: " << min_distance[i] << endl;
    }
    neighbors_file.close();

    return 0;
}
