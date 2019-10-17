#include "Library.h"

using namespace std;

class HashTable {
private:
    vector<vector<int>>* pTable;
    int Size;
public:
    HashTable(int TableSize);
    int Hash(int HashCode);
    void Insert(int HashCode, vector<int> Point);
    void Insert_to_Vertex(int Vertex, vector<int> Point);
    vector<vector<int>>* Search_Neighbors(int HashCode);
    vector<vector<int>>* Search_Neighbors_to_Vertex(int Vertex);
    ~HashTable();
};
