#include "lib/Library.h"
#include "lib/LSH_Functions.h"
#include "lib/Helper_Functions.h"
#include "lib/HashTable.h"

using namespace std;

int main(int argc, char* argv[]){
    if(argc != 3) {
      cout << "We need input_file AND query_file!" << endl;
      return -1;
    }
    /* variable declaration | k = 4 default value */
    int error_code, k = 4;
    /* vectors for the data and query points */
    vector<vector<int>> dataset;
    vector<vector<int>> searchset;

    /* read data set and query set and load them in vectors */
    error_code = Read_input_files(&dataset, &searchset, argv[1], argv[2]);
    if (error_code == -1) return -1;

    /* d-dimensional vectors */
    int d = dataset[0].size();
    /* compute window for all hash tables */
    int w = compute_window(dataset);

    /* do brute force to find actual NNs */
#ifdef BRUTE_FORCE
    brute_force(&dataset, &searchset);
#endif

    /* TODO: loop for L */
    vector<vector<int>> s;
    generate_shifts(&s, w, d, k);

    /* H of size (k, dataset.size()) */
    vector<vector<int>> hash_functions;
    /* projections of data */
    vector<vector<int>> a_projects;
    /* internal vector H for pushing */
    vector<int> H;
    /* loop for K */
    for (int i = 0; i < k; i++) {
        a_projects.clear();
        projections(&a_projects, dataset, &(s[i]), w, d);

        H.clear();
        compute_hash(&H, a_projects, d, k, w);
        hash_functions.push_back(H);
    } /* end for */
    /* compute the amplified hashes for every item */
    vector<string> amplified_g;
    amplify_hash(&amplified_g, &hash_functions, k);
    /* end for */

    /* Now that we have the hash codes, lets put them in the hash table */
    int TableSize = dataset.size()/16;
    HashTable* MyHashTable = new HashTable(TableSize);

    /* Insert all items inside the Hash Table */
    for (int i = 0; i < dataset.size(); i++){
        /* TODO: convert amplified_g can't fit in an int 32 bit */
        // MyHashTable->Insert(stoi(amplified_g[i]), dataset[i]);
    }

    /* TODO: do the same for queries, and put them inside the hash table */

    return 0;
}