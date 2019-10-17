#include "Library.h"
#include "Helper_Functions.h"

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

    m1 = dataset.size();
    m2 = dataset.size();
    double ** c = new double* [m1];
    for (int i = 0; i < m1; i++) {
        c[i] = new double [m2];
    }
    //DTW(&c, &dataset, &dataset);

    /* Free allocated space */
    for (int i = 0; i < m1; i++) {
        delete(c[i]);
    }
    delete c;

    return 0;
}