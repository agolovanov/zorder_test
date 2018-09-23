#include "morton.h"
#include <iostream>
using namespace std;

int main() {
    const int pow = 8;
    const int size = 1 << pow;

    auto encoder = morton_encoder(pow);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int morton = encoder.encode(i, j);
            cout << "(" << i << ", " << j << ") -> " << morton << endl;
        }
    }

    return 0;
}
