#include "Library.h"
#include "Helper_Functions.h"
#include "Traversals.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "We need input_file" << endl;
        return -1;
    }
    /* variable declaration | k = 4 default value */
    int k = 4, L = 5;
    int m1, m2, error_code;
    /* vectors for the data and query points */
    vector<vector<double*>> dataset;
    /* read data set and query set and load them in vectors */
    error_code = Read_curve_files(&dataset, argv[1]);
    if (error_code == -1) return -1;

    /* dataset sizes */
    m1 = dataset.size();
    m2 = dataset.size();
    /* allocate space */
    double ** c = new double* [m1];
    for (int i = 0; i < m1; i++) {
        c[i] = new double [m2];
    }

    /* compute Dynamic Time Warping */
    for(int i = 0; i < dataset.size(); i++) {
        for (int j = dataset.size() - 1; j > 0; j--) {
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