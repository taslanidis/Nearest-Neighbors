#include "LSH.h"
#include "LSH_Functions.h"
#include "BHC.h"
#include "BHC_Functions.h"
#include "Helper_Functions.h"
#include "HashTable.h"

using namespace std;

template void BHC <int>(vector<vector<int>>*, vector<vector<int>>*, int, int, int, int, int, double, vector<vector<int>>*, int**, double**, int**);
template void BHC <double>(vector<vector<double>>*, vector<vector<double>>*, int, int, int, int, double, double, vector<vector<int>>*, double**, double**, int**);

template <typename Point>
void BHC (vector<vector<Point>>* dataset, vector<vector<Point>>* searchset, int k, int dim, int M, int probes, Point w, double R, vector<vector<int>>* R_Neighbors, Point** min_distance, double** time, int** nearest_neighbor){
    int d_size = dataset->size();
    int s_size = searchset->size();
    /* d-dimensional vectors */
    int d = (*dataset)[0].size();
    /* Size of Hash Table */
    int TableSize = d_size / 8;
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
    /* map to assign every g function to 0 or 1 */
    map <int, int> dictionary;
    /* results */
    vector<vector<vector<vector<Point>>>> ANN;
    /* store neighbors from probes */
    vector<vector<vector<Point>>> Neighbors;
    /* vector for bonus to store neighbors with radius r */
    vector<int> Curr_R_Neighbors;

    cout << "\nComputing w ... " << endl;
    //w = 4*compute_window(dataset);
    cout << "Computed w : " << w << endl;

    /* computing big numbers */
    int capital_M = pow(2, 32/k);
    int m = 3;
    int * power = new int [d-1];
    for (int j = 0; j < d-1; j++)
        power[j] = moduloPow(m, j, capital_M);

    /* loop for L, to create L amplified functions g */
    for (int l = 0; l < dim; l++) {
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
        /* compute the amplified hashes for every data item */
        amplify_hash(&temp_g, &hash_functions, k);
        data_amplified_g.push_back(temp_g);

        /* clear hash functions for search set */
        vector<vector<int>>().swap(hash_functions);
        vector<int>().swap(temp_g);

        /* ------------------------ SEARCH SET ----------------------------*/
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
        /* compute the amplified hashes for every query item */
        amplify_hash(&temp_g, &hash_functions, k);
        query_amplified_g.push_back(temp_g);

        /* clear hash functions for search set */
        vector<vector<int>>().swap(hash_functions);
        vector<int>().swap(temp_g);
    }
    /* free memory from the exponents array */
    delete[] power;

    fill_dictionary(&dictionary, data_amplified_g);

    int hypercube_size = pow(2,dim);
    int vertex = 0;

    vector<vector<Point>> MyVerticesTable[hypercube_size];
    for (int i = 0; i < dataset->size(); i++) {
        vertex = calculate_vertex(data_amplified_g, dictionary, i);
        MyVerticesTable[vertex].push_back((*dataset)[i]);
    }

    int ham_dist = 0;
    for (int i = 0; i < searchset->size(); i++) {
        vertex = calculate_vertex(query_amplified_g, dictionary, i);
        Neighbors.push_back(MyVerticesTable[vertex]);
        for (int j = 0; j < hypercube_size && Neighbors.size() < probes; j++) {
            ham_dist = hammingDistance(vertex, j);
            if(ham_dist == 1){
                Neighbors.push_back(MyVerticesTable[j]);
            }
        }
        ANN.push_back(Neighbors);
        vector<vector<vector<Point>>>().swap(Neighbors);
    }

    int distance = 0;
    /* default metric L1 Manhattan */
    int Metric = 1;
    int computations = 0;
    for (int q = 0; q < searchset->size(); q++) {
        auto start = chrono::high_resolution_clock::now();
        /* find for every query its neighbor vertices */
        for (int n = 0; n < ANN[q].size(); n++) {
            computations = 0;
            for (int j = 0; j < ANN[q][n].size(); j++) {
                distance = dist(&ANN[q][n][j], &searchset->at(q), dataset->at(0).size(), Metric);
                if(distance <= R){
                    Curr_R_Neighbors.push_back(ANN[q][n][j][0]);
                }
                if (((distance < (*min_distance)[q]) || (*min_distance)[q] == -1) && computations < M) {
                    (*min_distance)[q] = distance;
                    (*nearest_neighbor)[q] = ANN[q][n][j][0]; //todo:used to have +1
                }
                computations++;
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
        vector<int>().swap(Curr_R_Neighbors);
    }
}
