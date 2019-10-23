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
    int k = 4, L = 5, d = 2, L_vec = 5;
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
    delta = 0.0005;


    /* compute window for all hash tables (try *4 or *10) */
    //double w = 4*compute_window(dataset);
    double w = 4*28; // computed w and its 28
    double max_element = 0.0;
    /* maximum number of polygonal curve points */
    int max_points = 0;

    /* do brute force to find actual NNs */
    cout << "Exhaustive Search using DTW. It might take a while ..." << endl;
    vector<double> TrueDistances;
    vector<double> TrueTimes;
    curves_brute_force(&dataset, &searchset, &TrueDistances, &TrueTimes);
    cout << "Found exact neighbors. Proceeding to hashing ..." << endl;

    /*  --------- TODO: Loop this L times and then dtw on those L nn sets -------- */
    /* ----------------------- HASHING with ORTHOGONAL GRID ---------------------- */
    /* orthogonal grid of size d */
    vector<double> orthogonal_grid;
    shift_grid(&orthogonal_grid, delta, d);

    /* vector for hashed curves */
    vector<vector<double*>> hashed_curves;
    vector<double*> temp_hash;
    vector<vector<double>> data_vectored_curves;
    vector<vector<double>> search_vectored_curves;
    vector<double> curve;

    /* ------------------ DATA SET hashing ----------------- */
    /* hash all dataset curves */
    for (int i = 0; i < dataset.size(); i++) {
        hash_curve(&temp_hash, &dataset[i], &orthogonal_grid, delta, d);
        hashed_curves.push_back(temp_hash);
        temp_hash.clear();
    }

    /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
     * the equivalent points in our new grid, that the polygonal curve projects */
    int elements = 0;
    max_points = 0;
    /* concat the 2d points in every h to make it from (x1,y1)(x2,y2) to x1,y1,x2,y2 */
    for (int i = 0; i < hashed_curves.size(); i++) {
        for (int j = 0; j < hashed_curves[i].size(); j++) {
            curve.push_back(hashed_curves[i][j][0]);
            curve.push_back(hashed_curves[i][j][1]);
            /* find max element */
            if (hashed_curves[i][j][0] > max_element) {
                max_element = hashed_curves[i][j][0];
            }
            if (hashed_curves[i][j][1] > max_element) {
                max_element = hashed_curves[i][j][1];
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
            curve.push_back(hashed_curves[i][j][0]);
            curve.push_back(hashed_curves[i][j][1]);
            if (hashed_curves[i][j][0] > max_element) {
                max_element = hashed_curves[i][j][0];
            }
            if (hashed_curves[i][j][1] > max_element) {
                max_element = hashed_curves[i][j][1];
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
            data_vectored_curves[i].push_back((double)2*max_element);
        }
    }
    for (int i = 0; i < search_vectored_curves.size(); i++) {
        while (search_vectored_curves[i].size() < max_points) {
            search_vectored_curves[i].push_back((double)2*max_element);
        }
    }
    /* ----------- end of padding ------------ */

    /* --- initializations for LSH --- */
    int distance = 0;
    double *min_distance = new double [searchset.size()];
    int *nearest_neighbor = new int [searchset.size()];
    double *time = new double [searchset.size()];
    double max_af = 0.0;
    double average_af = 0.0;
    double curr_fraction = 0.0;
    double average_time = 0.0;
    for (int i = 0; i < searchset.size(); i++) {
        min_distance[i] = INT_MAX;
        nearest_neighbor[i] = -1;
        time[i] = 0;
    }

   /*for (int i = 0; i < data_vectored_curves.size(); i++){
        for (int j = 0; j < data_vectored_curves[i].size(); j++){
            cout << data_vectored_curves[i][j] << " | ";
        }
        cout << endl;
        getchar();
    }*/
    /* ---------------- Hashing them again with LSH ------------------ */
    LSH(&data_vectored_curves, &search_vectored_curves, k, L_vec, w, &min_distance, &time, &nearest_neighbor);

    /* TODO: min distance vector is for the lsh hashed data, we will do DTW now on the real curves to find the
     * true distance between the approximate nearest neighbors found by lsh*/
    for (int i = 0; i < searchset.size(); i++) {
        /* -1 because nn starts from 1 to N */
        min_distance[i] = DTW(&searchset[i], &dataset[nearest_neighbor[i] - 1]);
    }

    /* Array for DTW metric */
    /* TODO: compare with DTW the L different neighbors for every q*/

    /* Results for every curve query */
    int computations = 0;
    for (int q = 0; q < searchset.size(); q++) {
        curr_fraction = (double) min_distance[q] / TrueDistances[q];
        if (curr_fraction > max_af) max_af = curr_fraction;
        average_af += curr_fraction;
        average_time += time[q];
    }

    /* --- RESULTS --- */
    average_af = average_af / searchset.size();
    average_time = average_time / searchset.size();
    cout << "Variables used: | k_vec = " << k << " | L = " << L << " | L_vec = " << L_vec << endl;
    cout << "MAX Approximation Fraction (LSH Distance / True Distance) = " << max_af << endl;
    cout << "Average Approximation Fraction (LSH Distance / True Distance) = " << average_af << endl;
    cout << "Average Time of LSH Distance Computation = " << average_time << endl;

    /* open file to write results */
    ofstream neighbors_file;
    neighbors_file.open ("./output/curves_grid_lsh.txt");
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