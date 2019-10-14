#include "lib/HashTable.h"

using namespace std;

HashTable::HashTable(int TableSize) {
    pTable = new vector<vector<int>> * [TableSize];
    Size = TableSize;
    for (int i = 0; i < TableSize; i++) {
        pTable[i] = NULL;
    }
}

int HashTable::Hash(int HashCode) {
    return HashCode % Size;
}

void HashTable::Insert(int HashCode, vector<int> Point){
    int hashed = Hash(HashCode);
    pTable[hashed]->push_back(Point);
}

vector<vector<int>>* HashTable::Search_Neighbors(int HashCode){
    int hashed = Hash(HashCode);
    return pTable[hashed];
}

HashTable::~HashTable(){
    /* TODO: free memory recursively */
    delete(pTable);
    return;
}
