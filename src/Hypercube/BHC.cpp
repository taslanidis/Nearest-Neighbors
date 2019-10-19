#include "Library.h"
#include "LSH_Functions.h"
#include "BHC_Functions.h"
#include "Helper_Functions.h"
#include "HashTable.h"

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
    error_code = Read_point_files(&dataset, &searchset, argv[1], argv[2]);
    if (error_code == -1) return -1;

    vector<int> TrueDistances;
    vector<double> TrueTimes;
    /* do brute force to find actual NNs */
//#ifdef BRUTE_FORCE
    brute_force(&dataset, &searchset, &TrueDistances, &TrueTimes);
//#endif

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
    /* amplified hash for dataset*/
    vector<vector<int>> data_amplified_g;
    /* amplified hash for searchset */
    vector<vector<int>> query_amplified_g;
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

    vector<vector<vector<int>>> ANNi;
    for (int i = 0; i < searchset.size(); i++) {
        vertex = calculate_vertex(query_amplified_g, dictionary, i);
//        cout << i << " : " << vertex << endl;
        ANNi.push_back(MyVerticesTable[vertex]);
    }
//    cout << ANNi.size() << endl << ANNi[0].size() << endl << ANNi[0][0].size() <<endl;

    int distance = 0;
    int *min_distance = new int [searchset.size()];
    int *nearest_neighbor = new int [searchset.size()];
    double *time = new double [searchset.size()];
    double max_af = 0.0;
    double average_af = 0.0;
    double curr_fraction = 0.0;
    double average_time = 0.0;

    /* initialize arrays */
    for (int i = 0; i < searchset.size(); i++) {
        min_distance[i] = INT_MAX;
        nearest_neighbor[i] = -1;
        time[i] = 0;
    }

//    if (query_amplified_g[i][q] == data_amplified_g[i][ANN[i][q][j][0]])
//    TODO: fix probes!!
    /* for every query */
    int calculations = 0;
    int total = 0;     //NO NEED - just checking..
    for (int q = 0; q < searchset.size(); q++) {
        /* for every hash table L */
        auto start = chrono::high_resolution_clock::now();
        /* for every vector in the same bucket (max M calculations) */
        calculations = 0;
        total = 0;
        for (int j = 0; j < ANNi[q].size() && calculations < M; j++) {
            if (check_compatibility(&query_amplified_g, &data_amplified_g, q, ANNi[q][j][0])) {
                distance = dist(&ANNi[q][j], &searchset[q], d);
                if (distance < min_distance[q]) {
                    min_distance[q] = distance;
                    nearest_neighbor[q] = ANNi[q][j][0] + 1;
                }
                calculations++;
            }
            total++;
        }
//        cout << "Calculations = " << calculations <<endl;
//        cout << "Total = " << total  << endl;
//        getchar();
        curr_fraction = (double) min_distance[q] / TrueDistances[q];
        if (curr_fraction > max_af) max_af = curr_fraction;
        average_af += curr_fraction;
        auto finish = chrono::high_resolution_clock::now();
        auto elapsed = finish - start;
        double time_elapsed = chrono::duration<double>(elapsed).count();
        average_time += time_elapsed;
        time[q] = time_elapsed;
    }
    average_af = average_af / s_size;
    average_time = average_time / s_size;
    cout << "Variables used: | k = " << k << " | dim =  " << dim << " | M = " << M << " | probes " << probes << " | w = " << w << " | " <<endl;
    cout << "MAX Approximation Fraction (LSH Distance / True Distance) = " << max_af << endl;
    cout << "Average Approximation Fraction (LSH Distance / True Distance) = " << average_af << endl;
    cout << "Average Time of LSH Distance Computation = " << setprecision(9) << showpoint << fixed << average_time << endl;

    /* print results */
    /* open file to write results */
    ofstream neighbors_file;
    neighbors_file.open ("./output/nneighbors_bhc.txt");
    for (int i = 0; i < searchset.size(); i++) {
        neighbors_file << "Query: " << i + 1 << endl;
        neighbors_file << "Nearest Neighbor: " << nearest_neighbor[i]<< endl;
        neighbors_file << "distanceLSH: " << min_distance[i] << endl;
        neighbors_file << "distanceTrue: " << TrueDistances[i] << endl;
        neighbors_file << "tLSH: " << setprecision(9) << showpoint << fixed << time[i] << endl;
        neighbors_file << "tTrue: " << setprecision(9) << showpoint << fixed << TrueTimes[i] << endl << endl;
    }
    neighbors_file.close();

    return 0;
}
