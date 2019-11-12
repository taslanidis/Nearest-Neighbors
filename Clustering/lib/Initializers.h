#include <string>

using namespace std;

class Initializer {
protected:
    int K;
public:
    Initializer(){}
    virtual void init() {}
    virtual string get_name() {}
};

class Random_Selection : public Initializer {
private:
    string name = "Random Selection";
public:
    Random_Selection(int K){this->K = K;}
    void init();
    string get_name();
};

class KMeans_plusplus : public Initializer {
private:
    string name = "K-Means++";
public:
    KMeans_plusplus(int K){this->K = K;}
    void init();
    string get_name();
};