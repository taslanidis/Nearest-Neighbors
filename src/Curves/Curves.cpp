#include "Library.h"
#include "Helper_Functions.h"
#include "Traversals.h"

using namespace std;

int main(int argc, char* argv[]) {
    /* check arguments */
    if (argc != 2) {
        cout << "We need input_file" << endl;
        return -1;
    }

    /* variable declaration | k = 4 default value */
    int k = 4, L = 5, d = 8;
    int m1, m2, error_code, min, grid_dim;
    /* vectors for the data and query points */
    vector<vector<double*>> dataset;
    /* read data set and query set and load them in vectors */
    error_code = Read_curve_files(&dataset, argv[1]);
    if (error_code == -1) return -1;

    /* dataset sizes */
    m1 = dataset.size();
    m2 = dataset.size();
    /* orthogonal d-dimensional grid */
    min = (m1 < m2) ? m1 : m2;
    grid_dim = 4*d*min - 1;

    /* -------- ORTHOGONAL GRID pre processing ------- */
    /* orthogonal grid of size d */
    vector<int> orthogonal_grid;
    /* create grid */
    create_grid(&orthogonal_grid, grid_dim);
    /* shift grid by random t*/
    shift_grid(&orthogonal_grid, grid_dim);

    /* allocate space */
    double ** c = new double* [m1];
    for (int i = 0; i < m1; i++) {
        c[i] = new double [m2];
    }

    /* compute Dynamic Time Warping */
    for(int i = 0; i < dataset.size(); i++) {
        for (int j = 100; j > 0; j--) {
            DTW(&c, &dataset[i], &dataset[j]);
        }
    }

    /* Free allocated space */
    for (int i = 0; i < m1; i++) {
        delete(c[i]);
    }
    delete c;

    return 0;
}