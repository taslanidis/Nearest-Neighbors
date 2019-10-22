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
    int k = 4, L = 5, d = 2, L_vec = 1, epsilon = 0.5;
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

    /* ------------------------ HASHING with RANDOM PROJECTIONS ----------------------- */

    /* Let max_points be the max length of all input and possible query curves */
    int max_points = 0;
    for (int i = 0; i < dataset.size(); i++) {
        if (2*dataset[i][0][1] > max_points) {
            max_points = 2*dataset[i][0][1];
        }
    }
    for (int i = 0; i < searchset.size(); i++) {
        if (2*dataset[i][0][1] > max_points) {
            max_points = 2*dataset[i][0][1];
        }
    }


    /* ------------------ RELEVANT traversals ----------------- */

    /* Let matrix G have real, independent, normally distributed elements ∼N(0,1)
     * and dimensions K × d, where K = −d*log(epsilon) / epsilon2 */
    vector<double> G;
    double K = -d*log(epsilon)/pow(epsilon,2);
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    /* uniformly values in [0,1) */
    normal_real_distribution<double> distribution (0, 1);
    for (int i = 0; i < K; i++) {
        G.push_back(distribution(generator));
    }

    /* We project the Ui’s to vector x = [G·U1|···|G·Uu] ∈ R^uK
     * by multiplying G with every point vector, then concatenating. */

    /* Find relevant traversals, all the pairs (Ui,Vi)
     * and we keep the indexes for every pair of curve */
    vector<vector<int*>> relevant_traversals;
    vector<int*> traversals;
    for (int i = 0; i < dataset.size(); i++) {
        for (int j = 0; j < dataset.size(); j++) {
            if (i != j) {
                Relevant_Traversals(&traversals, dataset[i].size(), dataset[j].size(), 0, 0);
                relevant_traversals.push_back(traversals);
                traversals.clear();
            }
        }
    }

    /* TODO: Multiply G with every point vector and then concatenate into a vector for lsh */


    /* ---------------- Hashing them again with LSH ------------------ */
    vector<vector<int>> data_amplified_g;
    vector<vector<int>> query_amplified_g;
    vector<vector<vector<vector<double>>>> ANN;
    //LSH(&data_vectored_curves, &search_vectored_curves, k, L_vec, &data_amplified_g, &query_amplified_g, &ANN);

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
