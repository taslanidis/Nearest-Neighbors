#include "Library.h"

using namespace std;

template <typename T>
T modpow(T base, T exp, T modulus) {
    base %= modulus;
    T result = 1;
    while (exp > 0) {
        if (exp & 1) result = (result * base) % modulus;
        base = (base * base) % modulus;
        exp >>= 1;
    }
    return result;
}

int compute_window(vector<vector<int>>);
void projections(vector<vector<int>>*, vector<vector<int>>*, vector<int>*, int, int);
void generate_shifts(vector<vector<int>>*, int, int, int);
void compute_hash(vector<int>*, vector<vector<int>>*, int, int, int);
void amplify_hash(vector<int>*, vector<vector<int>>*, int);

long long moduloMultiplication(long long a, long long b, long long mod);
int modulo (int a, int b);