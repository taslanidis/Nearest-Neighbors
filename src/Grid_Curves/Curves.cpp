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
    delta = 0.05;


    /* compute window for all hash tables (try *4 or *10) */
    //double w = 4*compute_window(dataset);
    double w = 4*620;

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

    /* maximum number of polygonal curve points */
    int max_points = 0;

    /* ------------------ DATA SET hashing ----------------- */
    /* hash all dataset curves */
    for (int i = 0; i < dataset.size(); i++) {
        if (2*dataset[i][0][1] > max_points) {
            max_points = 2*dataset[i][0][1];
        }
        hash_curve(&temp_hash, &dataset[i], &orthogonal_grid, delta, d);
        hashed_curves.push_back(temp_hash);
        temp_hash.clear();
    }

    /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
     * the equivalent points in our new grid, that the polygonal curve projects */

    /* concat the 2d points in every h to make it from (x1,y1)(x2,y2) to x1,y1,x2,y2 */
    for (int i = 0; i < hashed_curves.size(); i++) {
        for (int j = 0; j < hashed_curves[i].size(); j++) {
            curve.push_back(hashed_curves[i][j][0]);
            curve.push_back(hashed_curves[i][j][1]);
        }
        data_vectored_curves.push_back(curve);
        curve.clear();
        curve.shrink_to_fit();
    }

    /* pad special number > max coord */
    for (int i = 0; i < data_vectored_curves.size(); i++) {
        while (data_vectored_curves[i].size() < max_points) {
            data_vectored_curves[i].push_back(0.0);
        }
    }

    /* ------------------ SEARCH SET hashing ----------------- */
    max_points = 0;
    hashed_curves.clear();
    hashed_curves.shrink_to_fit();
    /* hash all curves */
    for (int i = 0; i < searchset.size(); i++) {
        if (2*searchset[i][0][1] > max_points) {
            max_points = 2*searchset[i][0][1];
        }
        hash_curve(&temp_hash, &searchset[i], &orthogonal_grid, delta, d);
        hashed_curves.push_back(temp_hash);
        temp_hash.clear();
    }

    /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
     * the equivalent points in our new grid, that the polygonal curve projects */

    /* concat the 2d points in every h to make it from (x1,y1)(x2,y2) to x1,y1,x2,y2 */
    for (int i = 0; i < hashed_curves.size(); i++) {
        for (int j = 0; j < hashed_curves[i].size(); j++) {
            curve.push_back(hashed_curves[i][j][0]);
            curve.push_back(hashed_curves[i][j][1]);
        }
        search_vectored_curves.push_back(curve);
        curve.clear();
        curve.shrink_to_fit();
    }

    /* pad special number > max coord */
    for (int i = 0; i < search_vectored_curves.size(); i++) {
        while (search_vectored_curves[i].size() < max_points) {
            search_vectored_curves[i].push_back(0.0);
        }
    }

    /* --- initializations for LSH --- */
    int distance = 0;
    int *min_distance = new int [searchset.size()];
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