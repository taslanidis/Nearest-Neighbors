#include "Library.h"

using namespace std;

class HashTable {
private:
    vector<vector<int>> ** pTable;
    int Size;
public:
    HashTable(int TableSize);
    int Hash(int HashCode);
    void Insert(int HashCode, vector<int> Point);
    vector<vector<int>>* Search_Neighbors(int HashCode);
    ~HashTable();
};
