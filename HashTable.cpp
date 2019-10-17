#include "lib/HashTable.h"

using namespace std;

HashTable::HashTable(int TableSize) {
    pTable = new vector<vector<int>> [TableSize];
    Size = TableSize;
}

int HashTable::Hash(int HashCode) {
    return HashCode % Size;
}

void HashTable::Insert(int HashCode, vector<int> Point){
    int hashed = Hash(HashCode);
    pTable[hashed].push_back(Point);
}

void HashTable::Insert_to_Vertex(int Vertex, vector<int> Point){
    pTable[Vertex].push_back(Point);
}

vector<vector<int>>* HashTable::Search_Neighbors(int HashCode){
    int hashed = Hash(HashCode);
    return &pTable[hashed];
}

vector<vector<int>>* HashTable::Search_Neighbors_to_Vertex(int Vertex){
    return &pTable[Vertex];
}

HashTable::~HashTable(){
    /* TODO: free memory recursively */
    delete(pTable);
    return;
}
