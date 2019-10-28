#include "LSH.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"
#include "HashTable.h"

using namespace std;

template void LSH <int>(vector<vector<int>>*, vector<vector<int>>*, int, int, int, double, vector<vector<int>>*, int**, double**, int**);
template void LSH <double>(vector<vector<double>>*, vector<vector<double>>*, int, int, double, double, vector<vector<int>>*, double**, double**, int**);

template <typename Point>
void LSH (vector<vector<Point>>* dataset, vector<vector<Point>>* searchset, int k, int L, Point w, double R, vector<vector<int>>* R_Neighbors, Point** min_distance, double** time, int** nearest_neighbor){
    int d_size = dataset->size();
    int s_size = searchset->size();
    /* d-dimensional vectors */
    int d = (*dataset)[0].size();
    /* Size of Hash Table */
    int TableSize = d_size / 8;
    HashTable <Point> **MyHashTable = new HashTable <Point>* [L];
    /* vector containing (k,d) shifts */
    vector<vector<double>> s;
    /* H of size (k, dataset.size()) */
    vector<vector<int>> hash_functions;
    /* projections of data */
    vector<vector<int>> a_projects;
    /* internal vector H for pushing */
    vector<int> H;
    /* amplified hash for dataset*/
    vector<vector<int>> data_amplified_g;
    /* amplified hash for search-set */
    vector<vector<int>> query_amplified_g;
    /* temporary g vector for push back */
    vector<int> temp_g;
    /* results */
    vector<vector<vector<vector<Point>>>> ANN;
    /* approximate nearest neighbors for each L iter */
    vector<vector<vector<Point>>> ANNi;
    /* vector for bonus */
    vector<int> Curr_R_Neighbors;

    /* computing big numbers */
    int M = pow(2, 32/k);
    int m = 3;
    int * power = new int [d-1];
    for (int j = 0; j < d-1; j++)
        power[j] = moduloPower(m, j, M);

    /* loop for L, to create L hash tables */
    for (int l = 0; l < L; l++) {
        /* generate the random shifts */
        generate_shifts(&s, w, d, k);

        /* ----------------------- DATA SET -------------------------------*/
        /* loop for K */
        for (int i = 0; i < k; i++) {
            projections(&a_projects, dataset, &(s[i]), w, d);

            compute_hash(&H, &a_projects, &power, d, k, w);
            hash_functions.push_back(H);

            /* clear only deletes data, but does not free the underlying storage
            so we use swap */
            vector<int>().swap(H);
            vector<vector<int>>().swap(a_projects);
        } /* end for */
        /* compute the amplified hashes for every item */
        amplify_hash(&temp_g, &hash_functions, k);
        data_amplified_g.push_back(temp_g);

        /* clear hash functions for search set */
        vector<vector<int>>().swap(hash_functions);
        vector<int>().swap(temp_g);

        /* Now that we have the hash codes, lets put them in the hash table
        *  Insert all items inside the Hash Table */
        MyHashTable[l] = new HashTable<Point>(TableSize);
        for (int i = 0; i < dataset->size(); i++) {
            MyHashTable[l]->Insert(data_amplified_g[l][i], (*dataset)[i]);
        }
        /* ------------------------ SEARCH SET ----------------------------*/
        /* do the same for the queries, and put them inside the hash table */
        /* loop for K */
        for (int i = 0; i < k; i++) {
            projections(&a_projects, searchset, &(s[i]), w, d);

            compute_hash(&H, &a_projects, &power, d, k, w);
            hash_functions.push_back(H);

            /* clear only deletes data, but does not free the underlying storage
            so we use swap */
            vector<int>().swap(H);
            vector<vector<int>>().swap(a_projects);
        } /* end for */
        /* compute the amplified hashes for every item */
        amplify_hash(&temp_g, &hash_functions, k);
        query_amplified_g.push_back(temp_g);

        /* clear hash functions for search set */
        vector<vector<int>>().swap(hash_functions);
        vector<int>().swap(temp_g);

        /* calculate approximate nearest neighbors */
        for (int i = 0; i < searchset->size(); i++) {
            //cout << "query g: " << query_amplified_g[l][i] << " | data g: " << data_amplified_g[l][i] << endl;
            ANNi.push_back(*MyHashTable[l]->Search_Neighbors(query_amplified_g[l][i]));
        }
        ANN.push_back(ANNi);
        /* clear hash functions and s for next iteration */
        vector<vector<vector<Point>>>().swap(ANNi);
        vector<vector<double>>().swap(s);
    }

    double distance;
    int Metric = 1; //default
    int computations = 0;
    for (int q = 0; q < searchset->size(); q++) {
        /* for every hash table L */
        auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < ANN.size(); i++) {
            /* for every vector in the same bucket (max 4*L calculations) */
            computations = 0;
            for (int j = 0; j < ANN[i][q].size(); j++) {
                if (computations == 25 && R == 0) break;
                //cout << query_amplified_g[i][q] << " | " << data_amplified_g[i][(int)ANN[i][q][j][0]] << endl;
                if (query_amplified_g[i][q] == data_amplified_g[i][(int)ANN[i][q][j][0]]) {
                    distance = dist(&ANN[i][q][j], &searchset->at(q), dataset->at(0).size(), Metric);
                    if(distance <= R){
                        Curr_R_Neighbors.push_back(ANN[i][q][j][0]);
                    }
                    if (((distance < (*min_distance)[q]) || (*min_distance)[q] == -1) && computations < 25) {
                        (*min_distance)[q] = distance;
                        (*nearest_neighbor)[q] = ANN[i][q][j][0];
                    }
                    computations++;
                }
            }
        }
        auto finish = chrono::high_resolution_clock::now();
        auto elapsed = finish - start;
        double time_elapsed = chrono::duration<double>(elapsed).count();
        (*time)[q] = time_elapsed;
        sort( Curr_R_Neighbors.begin(), Curr_R_Neighbors.end() );
        Curr_R_Neighbors.erase( unique( Curr_R_Neighbors.begin(), Curr_R_Neighbors.end() ), Curr_R_Neighbors.end() );
        vector<int>::iterator position = find(Curr_R_Neighbors.begin(), Curr_R_Neighbors.end(), (*nearest_neighbor)[q]);
        if (position != Curr_R_Neighbors.end())
            Curr_R_Neighbors.erase(position);
        R_Neighbors->push_back(Curr_R_Neighbors);

        /* clean underlying memory of the vector */
        vector<int>().swap(Curr_R_Neighbors);
    }

    /* clean memory */
    vector<vector<vector<vector<Point>>>>().swap(ANN);
    delete[] power;
    for (int l = 0; l < L; l++) {
        delete MyHashTable[l];
    }
    delete[] MyHashTable;

}
