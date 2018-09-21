#include "morton.h"
#include <iostream>
using namespace std;

int main() {
    const int size = 64;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int morton = encode(i, j, size);
            cout << "(" << i << ", " << j << ") -> " << morton << endl;
        }
    }

    return 0;
}
