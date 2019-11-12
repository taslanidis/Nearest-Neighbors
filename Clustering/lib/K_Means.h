#include "Initializers.h"

using namespace std;

template <class Point>
class K_Means {
private:
    Initializer* initializer;
    int K;
public:
    K_Means(int, string);
    ~K_Means();
};
