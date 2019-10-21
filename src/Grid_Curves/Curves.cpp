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
    int k = 4, L = 5, d = 2, L_vec = 1;
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

    /* ----------------------- HASHING with ORTHOGONAL GRID ---------------------- */
    /* orthogonal grid of size d */
    vector<int> orthogonal_grid;
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

    /* ---------------- Hashing them again with LSH ------------------ */
    /* TODO: now every h is a vector, and we will call lsh for those h */
    vector<vector<int>> data_amplified_g;
    vector<vector<int>> query_amplified_g;
    vector<vector<vector<vector<int>>>> ANN;
    LSH(&data_vectored_curves, &search_vectored_curves, k, L_vec, &data_amplified_g, &query_amplified_g, &ANN);

    /* TODO: store in table the lsh result */

    /* initializations */
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

    /* Array for DTW metric */
    cout << "Computed the hashes\nNow computing DTW. It might take a while ..." << endl;
    /* allocate space */
    double ** c = new double* [m1];
    for (int i = 0; i < m1; i++) {
        c[i] = new double [m2];
    }

    /* for every query */
    int computations = 0;
    for (int q = 0; q < searchset.size(); q++) {
        /* for every hash table L */
        auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < ANN.size(); i++) {
            /* for every vector in the same bucket (max 4*L calculations) */
            computations = 0;
            for (int j = 0; j < ANN[i][q].size() && computations < 4 * L; j++) {
                if (query_amplified_g[i][q] == data_amplified_g[i][(int)ANN[i][q][j][0]]) {
                    distance = DTW(&c, &dataset[(int)ANN[i][q][j][0]], &searchset[q]);
                    if (distance < min_distance[q]) {
                        min_distance[q] = distance;
                        nearest_neighbor[q] = ANN[i][q][j][0] + 1;
                    }
                    computations++;
                }
            }
        }
    }

    /* Free allocated space */
    for (int i = 0; i < m1; i++) {
        delete(c[i]);
    }
    delete c;

    return 0;
}