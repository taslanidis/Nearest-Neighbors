#include "LSH.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"
#include "HashTable.h"

using namespace std;

template void LSH <int>(vector<vector<int>>*, vector<vector<int>>*, int, int, int, double, vector<vector<int>>*, int**, double**, int**);
template void LSH <double>(vector<vector<double>>*, vector<vector<double>>*, int, int, double, double, vector<vector<int>>*, double**, double**, int**);

template <typename Point>
void LSH (vector<vector<Point>>* dataset, vector<vector<Point>>* searchset, int k, int L, Point w, double R, vector<vector<int>>* R_Neighbors, Point** min_distance, double** time, int** nearest_neighbor){

    /* Size of dataset and searchset */
    int d_size = dataset->size();
    int s_size = searchset->size();
    /* Size of Vectors (d-dimensional vectors) */
    int d = (*dataset)[0].size();
    /* Size of Hash Table */
    int TableSize = d_size / 8;
    HashTable <Point> **MyHashTable = new HashTable <Point>* [L];
    /* Vector containing shifts of size(k,d) */
    vector<vector<double>> s;
    /* Vector containing H of size (k, d_size) */
    vector<vector<int>> hash_functions;
    /* Vector containing projections of data */
    vector<vector<int>> a_projects;
    /* Internal vector H for pushing */
    vector<int> H;
    /* Amplified hash for dataset*/
    vector<vector<int>> data_amplified_g;
    /* Amplified hash for searchset */
    vector<vector<int>> query_amplified_g;
    /* Temporary g vector for push back */
    vector<int> temp_g;
    /* ANN results */
    vector<vector<vector<vector<Point>>>> ANN;
    /* ANN for each L iteration */
    vector<vector<vector<Point>>> ANNi;
    /* Vector for R-Neighbors (BONUS) */
    vector<int> Curr_R_Neighbors;

    /* Computing big numbers */
    int M = pow(2, 32/k);
    int m = 3;                              //chosen after testing because of the results that it gave us
    int * power = new int [d-1];
    for (int j = 0; j < d-1; j++)
        power[j] = moduloPower(m, j, M);

    /* Loop for creation of L hash tables */
    for (int l = 0; l < L; l++) {

        /* Shifts Generation */
        generate_shifts(&s, w, d, k);

        /* ----------------------- DATA SET -------------------------------*/

        /* Loop for creation of K hi to create amplified g*/
        for (int i = 0; i < k; i++) {
            /* Projections Computation */
            projections(&a_projects, dataset, &(s[i]), w, d);
            /* Hash Computation */
            compute_hash(&H, &a_projects, &power, d, k, w);
            hash_functions.push_back(H);

            /* clear only deletes data, but does not free the underlying storage
            so we use swap */
            vector<int>().swap(H);
            vector<vector<int>>().swap(a_projects);
        }

        /* Amplified hash computation */
        amplify_hash(&temp_g, &hash_functions, k);
        data_amplified_g.push_back(temp_g);

        vector<vector<int>>().swap(hash_functions);
        vector<int>().swap(temp_g);

        /* Insert all items inside the Hash Table */
        MyHashTable[l] = new HashTable<Point>(TableSize);
        for (int i = 0; i < dataset->size(); i++) {
            MyHashTable[l]->Insert(data_amplified_g[l][i], (*dataset)[i]);
        }

        /* ------------------------ SEARCH SET ----------------------------*/

        /* Loop for creation of K hi to create amplified g*/
        for (int i = 0; i < k; i++) {
            /* Projections Computation */
            projections(&a_projects, searchset, &(s[i]), w, d);
            /* Hash Computation */
            compute_hash(&H, &a_projects, &power, d, k, w);
            hash_functions.push_back(H);

            vector<int>().swap(H);
            vector<vector<int>>().swap(a_projects);
        }

        /* Amplified hash computation */
        amplify_hash(&temp_g, &hash_functions, k);
        query_amplified_g.push_back(temp_g);

        vector<vector<int>>().swap(hash_functions);
        vector<int>().swap(temp_g);

        /* ANN Calculation */
        for (int i = 0; i < searchset->size(); i++) {
            ANNi.push_back(*MyHashTable[l]->Search_Neighbors(query_amplified_g[l][i]));
        }
        ANN.push_back(ANNi);

        /* clear hash functions and s for next iteration */
        vector<vector<vector<Point>>>().swap(ANNi);
        vector<vector<double>>().swap(s);
    }

    double distance;
    /* default metric L1 Manhattan */
    int Metric = 1;
    int computations = 0;
    /* For every query in searchset */
    for (int q = 0; q < searchset->size(); q++) {
        auto start = chrono::high_resolution_clock::now();
        /* For every hash table L */
        for (int i = 0; i < ANN.size(); i++) {
            computations = 0;
            /* For every vector in the same bucket (max 4*L calculations - we use 25 as limit) */
            for (int j = 0; j < ANN[i][q].size(); j++) {
                if (computations == 25 && R == 0) break;
                /* Check points that have the same amplified g */
                if (query_amplified_g[i][q] == data_amplified_g[i][(int)ANN[i][q][j][0]]) {
                    distance = dist(&ANN[i][q][j], &searchset->at(q), dataset->at(0).size(), Metric);
                    /* Find R-Neighbors */
                    if(distance <= R){
                        Curr_R_Neighbors.push_back(ANN[i][q][j][0]);
                    }
                    /* Find Nearest Neighbor */
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
