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
    vector<vector<vector<vector<double>>>>** TraversalsTable = new vector<vector<vector<vector<double>>>>* [M];
    for (int i = 0; i < M; i++) {
        TraversalsTable[i] = new vector<vector<vector<vector<double>>>> [M];
    }

    vector<vector<vector<int>>> traversals;
    vector<vector<vector<vector<double>>>> traversals_to_coordinates;
    vector<vector<vector<double>>> traverse_with_coords;
    vector<vector<double>> pair_coords;
    vector<double> coords;
    int index_x;
    int index_y;
    /* for every pair of curves */
    for (int i = 0; i < dataset.size(); i++) {
        for (int j = 0; j < searchset.size(); j++) {
            /* find all relevant traversals */
            Relevant_Traversals(&traversals, dataset[i].size(), searchset[j].size());
            /* Convert the vector of pair with indexes, to a vector with coordinates */
            /* for all traversals of size */
            for (int k = 0; k < traversals.size(); k++) {
                /* for every pair */
                for (int t = 0; t < traversals[k].size(); k++) {
                    /* for U and V */
                    index_x = traversals[k][t][0];
                    index_y = traversals[k][t][1];
                    coords.push_back(dataset[i][index_x][0]);
                    coords.push_back(dataset[i][index_x][1]);
                    pair_coords.push_back(coords);
                    vector<double>().swap(coords);
                    coords.push_back(searchset[j][index_y][0]);
                    coords.push_back(searchset[j][index_y][1]);
                    pair_coords.push_back(coords);
                    traverse_with_coords.push_back(pair_coords);
                    vector<vector<double>>().swap(pair_coords);
                }
                traversals_to_coordinates.push_back(traverse_with_coords);
                vector<vector<vector<double>>>().swap(traverse_with_coords);
            }
            /* cell li,lj contains all relevant traversals of length li,lj curves. */
            for (int k = 0; k < traversals.size(); k++) {
                /* push back the traversal on the appropriate index (length of the curves) */
                TraversalsTable[(int)dataset[i][0][1]-1][(int)searchset[j][0][1]-1].push_back(traversals_to_coordinates[k]);
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
    vector<vector<double>> Vectored_Traversals_X;
    vector<vector<double>> Vectored_Traversals_Y;
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
                Vectored_Traversals_X.push_back(termU);
                Vectored_Traversals_Y.push_back(termV);
                vector<double>().swap(termU);
                vector<double>().swap(termV);
            }
        }
    }

    cout << "Computed the vectored traversals" << endl;


    /* ---------------- Hashing them again with LSH ------------------ */
    double w = 600;
    double R = 0;
    /* results */
    double *min_distance = new double[searchset.size()];
    int *nearest_neighbor = new int[searchset.size()];
    double *time = new double[searchset.size()];
    /* init arrays */
    for (int i = 0; i < searchset.size(); i++) {
        min_distance[i] = INT_MAX;
        nearest_neighbor[i] = -1;
        time[i] = 0;
    }
    /* results for bonus */
    vector<vector<int>> R_neighbors;

    cout << "Calling LSH ..." << endl;
    LSH(&Vectored_Traversals_X, &Vectored_Traversals_Y, k, L_vec, w, R, &R_neighbors, &min_distance, &time, &nearest_neighbor);

    cout << "LSH ended" << endl;

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

    return 0;
}
