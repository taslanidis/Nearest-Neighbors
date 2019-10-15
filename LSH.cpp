#include "lib/Library.h"
#include "lib/LSH_Functions.h"
#include "lib/Helper_Functions.h"
#include "lib/HashTable.h"

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 3) {
      cout << "We need input_file AND query_file!" << endl;
      return -1;
    }
    /* variable declaration | k = 4 default value */
    int error_code, k = 4, L = 5;
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
    /* Size of Hash Table */
    int TableSize = dataset.size() / 8;
    HashTable *MyHashTable[L];
    /* vector containing (k,d) shifts */
    vector<vector<int>> s;
    /* H of size (k, dataset.size()) */
    vector<vector<int>> hash_functions;
    /* projections of data */
    vector<vector<int>> a_projects;
    /* internal vector H for pushing */
    vector<int> H;
    /* amplified hash */
    vector <int> amplified_g;
    /* results */
    vector<vector<vector<vector<int>>>> ANN;
    /* loop for L, to create L hash tables */
    for (int l = 0; l < L; l++) {
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
        /* compute the amplified hashes for every item */
        amplify_hash(&amplified_g, &hash_functions, k);

        /* Now that we have the hash codes, lets put them in the hash table
        *  Insert all items inside the Hash Table */
        MyHashTable[l] = new HashTable(TableSize);
        for (int i = 0; i < dataset.size(); i++) {
            MyHashTable[l]->Insert(amplified_g[i], dataset[i]);
        }

        /* ------------------------ SEARCH SET ----------------------------*/
        /* do the same for the queries, and put them inside the hash table */
        /* loop for K */
        for (int i = 0; i < k; i++) {
            a_projects.clear();
            projections(&a_projects, searchset, &(s[i]), w, d);

            H.clear();
            compute_hash(&H, a_projects, d, k, w);
            hash_functions.push_back(H);
        } /* end for */
        /* compute the amplified hashes for every item */
        amplified_g.clear();
        amplify_hash(&amplified_g, &hash_functions, k);

        /* calculate approximate nearest neighbors */
        vector<vector<vector<int>>> ANNi;
        for (int i = 0; i < searchset.size(); i++) {
            ANNi.push_back(*MyHashTable[l]->Search_Neighbors(amplified_g[i]));
        }
        ANN.push_back(ANNi);
    }

    int distance = 0;
    int *min_distance = new int [searchset.size()];
    int *nearest_neighbor = new int [searchset.size()];

    /* initialize arrays */
    for (int i = 0; i < searchset.size(); i++) {
        min_distance[i] = INT_MAX;
        nearest_neighbor[i] = -1;
    }

    /* for every hash table L */
    for (int i = 0; i < ANN.size(); i++){
        /* for every query */
        for (int q = 0; q < searchset.size(); q++) {
            /* for every vector in the same bucket (max 3*L calculations) */
            for (int j = 0; j < ANN[i][q].size() && j < 4 * L; j++) {
                /* TODO: I have to check for same g(x) also */
                distance = dist(&ANN[i][q][j], &searchset[q], d);
                if (distance < min_distance[q]) {
                    min_distance[q] = distance;
                    nearest_neighbor[q] = ANN[i][q][j][0];
                }
            }
        }
    }

    /* print results */
    /* open file to write results */
    ofstream neighbors_file;
    neighbors_file.open ("nneighbors_lsh.txt");
    for (int i = 0; i < searchset.size(); i++) {
        neighbors_file << "Item: " << i + 1 << ", Neighbor: " << nearest_neighbor[i] << " | Distance: " << min_distance[i] << endl;
    }
    neighbors_file.close();

    return 0;
}