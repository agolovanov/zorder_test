/*
 * main.cpp
 *
 *  Created on: Sep 16, 2018
 *      Author: elrond16
 */

#include "array.h"
#include <iostream>
using namespace std;

int main(int argc, char **argv) {
    auto a = cached_array<double>(10, 10);

    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            cout << a(i, j) << " ";
        }
        cout << endl;
    }

    return 0;
}
