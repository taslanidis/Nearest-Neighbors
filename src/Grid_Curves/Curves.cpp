#include "Library.h"
#include "Helper_Functions.h"
#include "Traversals.h"
#include "LSH.h"

using namespace std;

int main(int argc, char* argv[]) {
    /* check arguments */
    if (argc != 3) {
        cout << "We need input_file and query file" << endl;
        return -1;
    }

    /* variable declaration | k = 4 default value */
    int k_vec = 2, L_grid = 4, d = 2, L_vec = 1;
    int m1, m2, error_code, min;
    double delta;
    /* vectors for the data and query points */
    vector<vector<double*>> dataset;
    vector<vector<double*>> searchset;
    /* read data set and query set and load them in vectors */
    error_code = Read_curve_files(&dataset, &searchset, argv[1], argv[2]);
    if (error_code == -1) return -1;

    /* dataset sizes */
    m1 = dataset.size();
    m2 = searchset.size();
    /* calculate delta */
    min = (m1 < m2) ? m1 : m2;
    //delta = 4*d*min - 1;
    delta = 0.15; // 0.15 gave us good results


    /* compute window for all hash tables (try *4 or *10) */
    //double w = 4*compute_window(dataset);
    double w = 10*60; // computed w and its ~28-32
    double max_element = 0.0;
    int max_points = 0, elements = 0;

    /* do brute force to find actual NNs */
    cout << "Exhaustive Search using DTW. It might take a while ..." << endl;
    vector<double> TrueDistances;
    vector<double> TrueTimes;
    //curves_brute_force(&dataset, &searchset, &TrueDistances, &TrueTimes);
    cout << "Found exact neighbors. Proceeding to hashing ..." << endl;

    /* orthogonal grid of size d */
    vector<double> orthogonal_grid;
    /* vector for hashed curves */
    vector <vector<double *>> hashed_curves;
    vector<double *> temp_hash;
    /* vectored curves */
    vector <vector<double>> data_vectored_curves;
    vector <vector<double>> search_vectored_curves;
    /* temp */
    vector<double> curve;
    /* RESULTS storing */
    double *min_distance;
    int *nearest_neighbor;
    double *time;
    vector<int*> lsh_neighbors;

    /*  ----------- Loop this L times and then dtw on those L nn sets -------- */
    for (int i = 0; i < L_grid; i++) {
        /* ----------------------- HASHING with ORTHOGONAL GRID ---------------------- */
        shift_grid(&orthogonal_grid, delta, d);

        /* ------------------ DATA SET hashing ----------------- */
        /* hash all dataset curves */
        for (int i = 0; i < dataset.size(); i++) {
            hash_curve(&temp_hash, &dataset[i], &orthogonal_grid, delta, d);
            hashed_curves.push_back(temp_hash);
            temp_hash.clear();
        }

        /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
         * the equivalent points in our new grid, that the polygonal curve projects */
        elements = 0;
        max_points = 0;
        /* concat the 2d points in every h to make it from (x1,y1)(x2,y2) to x1,y1,x2,y2 */
        for (int i = 0; i < hashed_curves.size(); i++) {
            for (int j = 0; j < hashed_curves[i].size(); j++) {
                curve.push_back(hashed_curves[i][j][0]);
                curve.push_back(hashed_curves[i][j][1]);
                /* find max element */
                /* j != 0 because at index 0 there is the id and length of curve */
                if (j != 0) {
                    if (hashed_curves[i][j][0] > max_element) {
                        max_element = hashed_curves[i][j][0];
                    }
                    if (hashed_curves[i][j][1] > max_element) {
                        max_element = hashed_curves[i][j][1];
                    }
                }
                elements++;
            }

            if (elements > max_points) {
                max_points = elements;
            }
            elements = 0;
            data_vectored_curves.push_back(curve);
            curve.clear();
            curve.shrink_to_fit();
        }

        /* ------------------ SEARCH SET hashing ----------------- */
        hashed_curves.clear();
        hashed_curves.shrink_to_fit();
        /* hash all curves */
        for (int i = 0; i < searchset.size(); i++) {
            hash_curve(&temp_hash, &searchset[i], &orthogonal_grid, delta, d);
            hashed_curves.push_back(temp_hash);
            temp_hash.clear();
        }

        /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
         * the equivalent points in our new grid, that the polygonal curve projects */

        /* concat the 2d points in every h to make it from (x1,y1)(x2,y2) to x1,y1,x2,y2 */
        elements = 0;
        for (int i = 0; i < hashed_curves.size(); i++) {
            for (int j = 0; j < hashed_curves[i].size(); j++) {
                /* we always push back the first element - even the id on (0,0)
                 * we won't push back the dimensions of a curve in (0,1) */
                curve.push_back(hashed_curves[i][j][0]);
                /* j != 0 because at index 0 there is the id and length of curve */
                if (j != 0) {
                    /* we dont push back the length size for the lsh */
                    curve.push_back(hashed_curves[i][j][1]);
                    if (hashed_curves[i][j][0] > max_element) {
                        max_element = hashed_curves[i][j][0];
                    }
                    if (hashed_curves[i][j][1] > max_element) {
                        max_element = hashed_curves[i][j][1];
                    }
                }
                elements++;
            }
            if (elements > max_points) {
                max_points = elements;
            }
            elements = 0;
            search_vectored_curves.push_back(curve);
            curve.clear();
            curve.shrink_to_fit();
        }

        /* ----------------- PADDING for both sets ------------ */
        /* pad special number > max coord */
        for (int i = 0; i < data_vectored_curves.size(); i++) {
            while (data_vectored_curves[i].size() < max_points) {
                data_vectored_curves[i].push_back(2 * max_element);
            }
        }
        for (int i = 0; i < search_vectored_curves.size(); i++) {
            while (search_vectored_curves[i].size() < max_points) {
                search_vectored_curves[i].push_back(2 * max_element);
            }
        }
        /* ----------- end of padding ------------ */

        /* --- initializations for LSH --- */
        min_distance = new double[searchset.size()];
        nearest_neighbor = new int[searchset.size()];
        time = new double[searchset.size()];

        for (int i = 0; i < searchset.size(); i++) {
            min_distance[i] = -1;
            nearest_neighbor[i] = -1;
            time[i] = 0;
        }


         /*for (int i = 0; i < search_vectored_curves.size(); i++){
             for (int j = 0; j < search_vectored_curves[i].size(); j++){
                 cout << search_vectored_curves[i][j] << " | ";
             }
             cout << endl;
             getchar();
         }*/

        /* ---------------- Hashing them again with LSH ------------------ */
        LSH(&data_vectored_curves, &search_vectored_curves, k_vec, L_vec, w, &min_distance, &time, &nearest_neighbor);

        /* store results for all iterations of hashing */
        lsh_neighbors.push_back(nearest_neighbor);

        /* clean vectors for next iteration */
        orthogonal_grid.clear();
        hashed_curves.clear();
        temp_hash.clear();
        data_vectored_curves.clear();
        search_vectored_curves.clear();
        curve.clear();
    }
    cout << "lsh ended" << endl;

    /* min distance vector is for the lsh hashed data, we will do DTW now on the real curves to find the
     * true distance between the approximate nearest neighbors found by lsh*/
    double distance = 0.0;
    double max_af = 0.0;
    double average_af = 0.0;
    double curr_fraction = 0.0;
    double average_time = 0.0;

    min_distance = new double[searchset.size()];
    nearest_neighbor = new int[searchset.size()];
    for (int i = 0; i < searchset.size(); i++) {
        min_distance[i] = INT_MAX;
        nearest_neighbor[i] = -1;
    }


    for (int i = 0; i < searchset.size(); i++) {
       for (int j = 0; j < L_grid; j++) {
           if (lsh_neighbors[j][i] == -1) continue;
           distance = DTW(&searchset[i], &dataset[lsh_neighbors[j][i]]);
           if (j == 0) {
               min_distance[i] = distance;
               nearest_neighbor[i] = lsh_neighbors[j][i];
           } else if (distance < min_distance[i]) {
               min_distance[i] = distance;
               nearest_neighbor[i] = lsh_neighbors[j][i];
           }
       }
    }

//    /* compare with DTW the L different neighbors for every q*/
//    /* Results for every curve query */
//    int computations = 0;
//    for (int q = 0; q < searchset.size(); q++) {
//        curr_fraction = (double) min_distance[q] / TrueDistances[q];
//        if (curr_fraction > max_af) max_af = curr_fraction;
//        average_af += curr_fraction;
//    }
//
//    /* --- RESULTS --- */
//    average_af = average_af / searchset.size();
//    average_time = average_time / searchset.size();
//    cout << "Variables used: | k_vec = " << k_vec << " | L_grid = " << L_grid << " | L_vec = " << L_vec << endl;
//    cout << "MAX Approximation Fraction (LSH Distance / True Distance) = " << max_af << endl;
//    cout << "Average Approximation Fraction (LSH Distance / True Distance) = " << average_af << endl;

    /* open file to write results */
    ofstream neighbors_file;
    neighbors_file.open ("./output/curves_grid_lsh.txt");
    for (int i = 0; i < searchset.size(); i++) {
        neighbors_file << "Query: " << i + 1 << endl;
        neighbors_file << "Method: LSH" << endl;
        neighbors_file << "HashFunction: LSH" << endl;
        if (nearest_neighbor[i] != -1) {
            neighbors_file << "Found Nearest Neighbor: " << nearest_neighbor[i] << endl;
            neighbors_file << "distanceFound: " << min_distance[i] << endl << endl;
        }else{
            neighbors_file << "Found Nearest Neighbor: Fail" << endl << endl;
        }
    }
    neighbors_file.close();

    return 0;
}