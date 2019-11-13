#include <string>
#include <vector>

using namespace std;

template <class Point>
class Assigner {
public:
    Assigner(){}
    virtual void assign(vector<vector<Point>>*) {}
    virtual string get_name() {}
};

template <class Point>
class Lloyd_assignment : public Assigner<Point> {
private:
    string name = "Lloyd's Assignment";
public:
    Lloyd_assignment(){};
    void assign(vector<vector<Point>>*);
    string get_name();
};

template <class Point>
class Inverse_assignment : public Assigner<Point> {
private:
    string name = "Inverse Assignment";
public:
    Inverse_assignment(){};
    void assign(vector<vector<Point>>*);
    string get_name();
};