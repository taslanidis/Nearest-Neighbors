#include "LSH.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"
#include "HashTable.h"

using namespace std;

/* TODO: LSH template types */

void LSH (vector<vector<int>>* dataset, vector<vector<int>>* searchset, int k, int L, vector<vector<int>>* data_amplified_g, vector<vector<int>>* query_amplified_g, vector<vector<vector<vector<int>>>>* ANN){
    int d_size = dataset->size();
    int s_size = searchset->size();
    /* d-dimensional vectors */
    int d = (*dataset)[0].size();
    /* compute window for all hash tables (try *4 or *10) */
    //int w = 4*compute_window(dataset);
    int w = 4*1164;
    /* Size of Hash Table */
    int TableSize = d_size / 8;
    HashTable *MyHashTable[L];
    /* vector containing (k,d) shifts */
    vector<vector<int>> s;
    /* H of size (k, dataset.size()) */
    vector<vector<int>> hash_functions;
    /* projections of data */
    vector<vector<int>> a_projects;
    /* internal vector H for pushing */
    vector<int> H;
    /* temporary g vector for push back */
    vector<int> temp_g;
    /* loop for L, to create L hash tables */
    for (int l = 0; l < L; l++) {
        /* generate the random shifts */
        generate_shifts(&s, w, d, k);

        /* ----------------------- DATA SET -------------------------------*/
        /* loop for K */
        for (int i = 0; i < k; i++) {
            projections(&a_projects, dataset, &(s[i]), w, d);

            compute_hash(&H, &a_projects, d, k, w);
            hash_functions.push_back(H);
            H.clear();
            H.shrink_to_fit();
            a_projects.clear();
            a_projects.shrink_to_fit();
        } /* end for */
        /* compute the amplified hashes for every item */
        amplify_hash(&temp_g, &hash_functions, k);
        data_amplified_g->push_back(temp_g);

        /* Now that we have the hash codes, lets put them in the hash table
        *  Insert all items inside the Hash Table */
        MyHashTable[l] = new HashTable(TableSize);
        for (int i = 0; i < dataset->size(); i++) {
            MyHashTable[l]->Insert(temp_g[i], (*dataset)[i]);
        }

        /* clear hash functions for search set */
        hash_functions.clear();
        hash_functions.shrink_to_fit();
        temp_g.clear();
        temp_g.shrink_to_fit();

        /* ------------------------ SEARCH SET ----------------------------*/
        /* do the same for the queries, and put them inside the hash table */
        /* loop for K */
        for (int i = 0; i < k; i++) {
            projections(&a_projects, searchset, &(s[i]), w, d);

            compute_hash(&H, &a_projects, d, k, w);
            hash_functions.push_back(H);
            H.clear();
            H.shrink_to_fit();
            a_projects.clear();
            a_projects.shrink_to_fit();
        } /* end for */
        /* compute the amplified hashes for every item */
        amplify_hash(&temp_g, &hash_functions, k);
        query_amplified_g->push_back(temp_g);

        /* calculate approximate nearest neighbors */
        vector<vector<vector<int>>> ANNi;
        for (int i = 0; i < searchset->size(); i++) {
            ANNi.push_back(*MyHashTable[l]->Search_Neighbors(temp_g[i]));
        }
        ANN->push_back(ANNi);

        /* clear hash functions and s for next iteration */
        hash_functions.clear();
        hash_functions.shrink_to_fit();
        temp_g.clear();
        temp_g.shrink_to_fit();
        s.clear();
        s.shrink_to_fit();
    }
}