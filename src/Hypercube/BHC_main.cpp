#include "Library.h"
#include "LSH.h"
#include "LSH_Functions.h"
#include "BHC.h"
#include "BHC_Functions.h"
#include "Helper_Functions.h"
#include "HashTable.h"

using namespace std;

int main(int argc, char* argv[]) {
    /* Arguments Check */
    int k = 4, dim = 3, M = 10, probes = 2;              //Default values: k is the number of hi concatenated to form g - dim is number of hypercube's vertices
    string data_file, search_file, results_file;
    int df = 0, sf = 0, rf = 0;
    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if(arg == "-d"){
            if(i+1 < argc){
                data_file = argv[++i];
                if(access( data_file.c_str(), F_OK ) == -1){
                    cout << "-- Wrong <input file> -- \n";
                    show_bhc_usage(argv[0]);
                    return -1;
                }
                df = 1;
            }else{
                cout << "-- NO <input file> -- \n";
                show_bhc_usage(argv[0]);
                return -1;
            }
        }else if (arg == "-q"){
            if(i+1 <= argc){
                search_file = argv[++i];
                if(access( search_file.c_str(), F_OK ) == -1){
                    cout << "-- Wrong <query file> -- \n";
                    show_bhc_usage(argv[0]);
                    return -1;
                }
                sf = 1;
            }else{
                cout << "-- NO <query file> -- \n";
                show_bhc_usage(argv[0]);
                return -1;
            }
        }else if(arg == "-o"){
            if(i+1 <= argc){
                results_file = argv[++i];
                rf = 1;
            }else{
                cout << "-- NO <output file> -- \n";
                show_bhc_usage(argv[0]);
                return -1;
            }
        }else if(arg == "-k"){
            dim = atoi(argv[++i]);
        }else if(arg == "-M") {
            M = atoi(argv[++i]);
        }else if(arg == "-probes"){
            probes = atoi(argv[++i]);
        }else{
            show_bhc_usage(argv[0]);
            return -1;
        }
    }
    if(df == 0 || sf == 0 || rf == 0){
        cout << "-- WRONG NUMBER OF FILES -- \n";
        show_bhc_usage(argv[0]);
        return -1;
    }

    /* vectors for the data and query points */
    vector<vector<int>> dataset;
    vector<vector<int>> searchset;

    /* read data set and query set and load them in vectors */
    int error_code = Read_point_files(&dataset, &searchset, data_file, search_file);
    if (error_code == -1) return -1;

    /* compute window for all hash tables (try *4 or *10) */
    //Point w = 4*compute_window(dataset);
    int w = 4*1140;

    /* do brute force to find actual NNs */
    vector<int> TrueDistances;
    vector<double> TrueTimes;
    // TODO: make compilation of brute force in makefile
    brute_force(&dataset, &searchset, &TrueDistances, &TrueTimes);

    /* results */
    int *min_distance = new int [searchset.size()];
    int *nearest_neighbor = new int [searchset.size()];
    double *time = new double [searchset.size()];
    for (int i = 0; i < searchset.size(); i++) {
        min_distance[i] = INT_MAX;
        nearest_neighbor[i] = -1;
        time[i] = 0;
    }

    /* ---- CALL BHC ---- */
    BHC(&dataset, &searchset, k, dim, M, probes, w,  &min_distance, &time, &nearest_neighbor);

    /* variables */
    double max_af = 0.0;
    double average_af = 0.0;
    double curr_fraction = 0.0;
    double average_time = 0.0;

    /* Results for every curve query */
    for (int q = 0; q < searchset.size(); q++) {
        curr_fraction = (double) min_distance[q] / TrueDistances[q];
        if (curr_fraction > max_af) max_af = curr_fraction;
        average_af += curr_fraction;
        average_time += time[q];
    }

    /* print results */
    average_af = average_af / searchset.size();
    average_time = average_time / searchset.size();
    cout << "Variables used: | k = " << k << " | dim = " << dim << " | M = " << M << " | probes = " << probes << endl;
    cout << "MAX Approximation Fraction (LSH Distance / True Distance) = " << max_af << endl;
    cout << "Average Approximation Fraction (LSH Distance / True Distance) = " << average_af << endl;
    cout << "Average Time of LSH Distance Computation = " << setprecision(9) << showpoint << fixed << average_time << endl;

    /* open file to write results */
    ofstream neighbors_file;
    neighbors_file.open ("./output/"+results_file);
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
