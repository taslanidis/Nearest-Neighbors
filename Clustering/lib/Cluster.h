#include "Initializers.h"

using namespace std;

template <class Point>
class Cluster {
private:
    Initializer<Point>* initializer;
    int K;
public:
    Cluster(int, string);
    void fit(vector<vector<Point>>*);
    ~Cluster();
};
