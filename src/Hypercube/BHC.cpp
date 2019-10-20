#include "LSH.h"
#include "LSH_Functions.h"
#include "BHC.h"
#include "BHC_Functions.h"
#include "Helper_Functions.h"
#include "HashTable.h"

using namespace std;

void BHC (vector<vector<int>> dataset, vector<vector<int>> searchset, int k, int dim, int M, int probes, vector<vector<int>> data_amplified_g, vector<vector<int>> query_amplified_g, vector<vector<vector<vector<int>>>>* ANN){
    int d_size = dataset.size();
    int s_size = searchset.size();
    /* d-dimensional vectors */
    int d = dataset[0].size();
    /* compute window for all hash tables (try *4 or *10) */
    //int w = 4*compute_window(dataset);
    int w = 4*1164;
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
    /* map to assign every g function to 0 or 1 */
    map <int, int> dictionary;
    /* loop for L, to create L amplified functions g */
    for (int l = 0; l < dim; l++) {
    /* generate the random shifts */
    generate_shifts(&s, w, d, k);

    /* ----------------------- DATA SET -------------------------------*/
    /* loop for K */
    for (int i = 0; i < k; i++) {
        projections(&a_projects, &dataset, &(s[i]), w, d);

        compute_hash(&H, &a_projects, d, k, w);
        hash_functions.push_back(H);
        H.clear();
        H.shrink_to_fit();
        a_projects.clear();
        a_projects.shrink_to_fit();
        } /* end for */
        /* compute the amplified hashes for every data item */
        amplify_hash(&temp_g, &hash_functions, k);
        data_amplified_g.push_back(temp_g);

        temp_g.clear();
        temp_g.shrink_to_fit();
        hash_functions.clear();
        hash_functions.shrink_to_fit();

        /* ------------------------ SEARCH SET ----------------------------*/
        /* loop for K */
        for (int i = 0; i < k; i++) {
            projections(&a_projects, &searchset, &(s[i]), w, d);

            compute_hash(&H, &a_projects, d, k, w);
            hash_functions.push_back(H);
            H.clear();
            H.shrink_to_fit();
            a_projects.clear();
            a_projects.shrink_to_fit();
        } /* end for */
        /* compute the amplified hashes for every query item */
        amplify_hash(&temp_g, &hash_functions, k);
        query_amplified_g.push_back(temp_g);

        temp_g.clear();
        temp_g.shrink_to_fit();
        hash_functions.clear();
        hash_functions.shrink_to_fit();
    }

    fill_dictionary(&dictionary, data_amplified_g);
    //    for(auto& x: dictionary){
    //        cout << x.first << ", " << x.second << endl;
    //    }

    int hypercube_size = pow(2,dim);
    int vertex = 0;

    vector<vector<int>> MyVerticesTable[hypercube_size];
    for (int i = 0; i < dataset.size(); i++) {
        vertex = calculate_vertex(data_amplified_g, dictionary, i);
        //        cout << i << " : " << vertex << endl;
        MyVerticesTable[vertex].push_back(dataset[i]);
    }

    vector<vector<vector<int>>> Neighbors;
    int ham_dist = 0;
    for (int i = 0; i < searchset.size(); i++) {
        vertex = calculate_vertex(query_amplified_g, dictionary, i);
        //        cout << i << " : " << vertex << endl;
        Neighbors.push_back(MyVerticesTable[vertex]);
        for (int j = 0; j < hypercube_size && Neighbors.size() < probes; j++) {
            ham_dist = hammingDistance(vertex, j);
            if(ham_dist == 1){
                Neighbors.push_back(MyVerticesTable[j]);
            }
        }
        ANN->push_back(Neighbors);
        Neighbors.clear();
        Neighbors.shrink_to_fit();
    }
}