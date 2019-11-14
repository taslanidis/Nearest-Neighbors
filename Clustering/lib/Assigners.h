#include <string>
#include <vector>
#include "Library.h"

using namespace std;

template <class Point>
class Assigner {
public:
    Assigner(){}
    virtual vector<vector<int>>* assign(vector<vector<Point>>*, vector<int>*) {return NULL;}
    virtual string get_name() {}
};

template <class Point>
class Lloyd_assignment : public Assigner<Point> {
private:
    string name = "Lloyd's Assignment";
public:
    Lloyd_assignment(){};
    vector<vector<int>>* assign(vector<vector<Point>>*, vector<int>*);
    string get_name();
};

template <class Point>
class Inverse_assignment : public Assigner<Point> {
private:
    string name = "Inverse Assignment";
public:
    Inverse_assignment(){};
    vector<vector<int>>* assign(vector<vector<Point>>*, vector<int>*);
    string get_name();
};