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
    double epsilon = 0.5;
    int m1, m2, error_code, min;
    double delta;
    /* vectors for the data and query points */
    vector<vector<double*>> dataset;
    vector<vector<double*>> searchset;
    /* read data set and query set and load them in vectors */
    error_code = Read_curve_files_max_dim(&dataset, &searchset, argv[1], argv[2], 15.00);
    if (error_code == -1) return -1;

    /* dataset sizes */
    m1 = dataset.size();
    m2 = searchset.size();
    /* calculate delta */
    min = (m1 < m2) ? m1 : m2;
    //delta = 4*d*min - 1;
    delta = 0.15;

    /* ------------------------ HASHING with RANDOM PROJECTIONS ----------------------- */

    /* Let M be the max length of all input and possible query curves */
    int M = 0;
    for (int i = 0; i < dataset.size(); i++) {
        if (dataset[i][0][1] > M) {
            M = dataset[i][0][1];
        }
    }
    for (int i = 0; i < searchset.size(); i++) {
        if (searchset[i][0][1] > M) {
            M = searchset[i][0][1];
        }
    }


    /* ------------------ RELEVANT traversals ----------------- */

    cout << "Creating MxM array where M is: " << M << endl;
    /* Find relevant traversals, all the pairs (Ui,Vi)
     * and we keep the indexes for every pair of curve */
    vector<vector<vector<int>>>** TraversalsTable = new vector<vector<vector<int>>>* [M];
    for (int i = 0; i < M; i++) {
        TraversalsTable[i] = new vector<vector<vector<int>>> [M];
    }

    vector<vector<vector<int>>> traversals;
    /* for every pair of curves */
    for (int i = 0; i < dataset.size(); i++) {
        for (int j = 0; j < searchset.size(); j++) {
            /* find all relevant traversals */
            Relevant_Traversals(&traversals, dataset[i].size(), searchset[j].size());
            /* cell li,lj contains all relevant traversals of length li,lj curves. */
            for (int k = 0; k < traversals.size(); k++) {
                /* push back the traversal on the appropriate index (length of the curves) */
                TraversalsTable[(int)dataset[i][0][1]-1][(int)searchset[j][0][1]-1].push_back(traversals[k]);
            }
        }
    }

    cout << "Array MxM is ready" << endl;

    /* Let matrix G have real, independent, normally distributed elements ∼N(0,1)
     * and dimensions K × d, where K = −d*log(epsilon) / epsilon2 */  //TODO: G must be Kxd
    double K = -d*log2(epsilon)/pow(epsilon,2);
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    /* uniformly values in [0,1) */
    vector<vector<double>> G;
    vector<double> G_temp;
    normal_distribution<double> distribution (0, 1);
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < d; j++) {
            G_temp.push_back(distribution(generator));
        }
        G.push_back(G_temp);
        G_temp.clear();
        G_temp.shrink_to_fit();
    }

    cout << "Vector G is ready" << endl;

    /* We project the Ui’s to vector x = [G·U1|···|G·Uu] ∈ R^uK
     * by multiplying G with every point vector, then concatenating. */
    vector<double> termU;
    vector<double> termV;
    vector<vector<double>> dataset;
    vector<vector<double>> searchset;
    /* for 1st dim of MxM table */
    for (int i = 0; i < M; i++) {
        /* for 2nd dim of MxM  table */
        for (int j = 0; j < M; j++) {
            /* For all traversals of size (i,j) */
            for (int h = 0; h < TraversalsTable[i][j].size(); h++) {
                /* For all pairs of traversal h */
                for (int t = 0; t < TraversalsTable[i][j][h].size(); t++) {
                    vector<double> U = TraversalsTable[i][j][h][t][0];
                    vector<double> V = TraversalsTable[i][j][h][t][1];
                    for (int k = 0; k < G.size(); k++){
                        for (int dim = 0; dim < d; dim++) {
                            termU.push_back(G[k][d] * U[dim]);
                            termV.push_back(G[k][d] * V[dim]);
                        }
                    }
                }
                dataset.push_back(termU);
                searchset.push_back(termV);
                vector<double>().swap(termU);
                vector<double>().swap(termV);
            }
        }
    }

    /* TODO: create M × M table: cell i, j contains all relevant traversals of length i, j curves */
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
                    distance = DTW(&dataset[(int)ANN[i][q][j][0]], &searchset[q]);          //todo: check this again,  there was a &c before
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
