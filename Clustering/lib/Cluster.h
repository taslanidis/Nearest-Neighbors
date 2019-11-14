#include "Initializers.h"
#include "Assigners.h"

using namespace std;

template <class Point>
class Cluster {
private:
    Initializer<Point>* initializer;
    Assigner<Point>* assigner;
    int K;
    vector<int>* centroids;
    vector<vector<int>>* clusters;
public:
    Cluster(int, string, string);
    void fit(vector<vector<Point>>*);
    ~Cluster();
};
