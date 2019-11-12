#include <string>
#include <vector>

using namespace std;

class Initializer {
protected:
    int K;
public:
    Initializer(){}
    virtual void init() {}
    virtual string get_name() {}
};

template <class Point>
class Random_Selection : public Initializer {
private:
    string name = "Random Selection";
public:
    Random_Selection(int K){this->K = K;}
    void init(vector<vector<Point>>*);
    string get_name();
};

template <class Point>
class KMeans_plusplus : public Initializer {
private:
    string name = "K-Means++";
public:
    KMeans_plusplus(int K){this->K = K;}
    void init(vector<vector<Point>>*);
    string get_name();
};