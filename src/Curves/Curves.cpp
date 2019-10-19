#include "Library.h"
#include "Helper_Functions.h"
#include "Traversals.h"

using namespace std;

int main(int argc, char* argv[]) {
    /* check arguments */
    if (argc != 3) {
        cout << "We need input_file and query file" << endl;
        return -1;
    }

    /* variable declaration | k = 4 default value */
    int k = 4, L = 5, d = 2;
    int m1, m2, error_code, min, delta;
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
    delta = 1;

    /* ------------ HASHING with ORTHOGONAL GRID ----------- */
    /* orthogonal grid of size d */
    vector<int> orthogonal_grid;
    /* create and shift grid by random vector t */
    vector<double> t;
    shift_grid(&orthogonal_grid, delta, d);

    /* vector for hashed curves */
    vector<vector<int*>> hashed_curves;
    vector<int*> temp_hash;
    /* hash all curves */
    for (int i = 0; i < dataset.size(); i++) {
        hash_curve(&temp_hash, &dataset[i], &orthogonal_grid, delta, d);
        hashed_curves.push_back(temp_hash);
        temp_hash.clear();
    }
    /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
     * the equivalent points in our new grid, that the polygonal curve projects */

    /* TODO: smth with the hashes to find the new points by adding the orth grid */

    /* ----------------- DTW ----------------- */
    cout << "Computed the hashes\nNow computing DTW. It might take a while ..." << endl;
    /* allocate space */
    double ** c = new double* [m1];
    for (int i = 0; i < m1; i++) {
        c[i] = new double [m2];
    }

    /* compute Dynamic Time Warping */
    for(int i = 0; i < dataset.size(); i++) {
        for (int j = 0; j < searchset.size(); j++) {
            DTW(&c, &dataset[i], &searchset[j]);
        }
    }

    /* Free allocated space */
    for (int i = 0; i < m1; i++) {
        delete(c[i]);
    }
    delete c;

    return 0;
}