/*
 * main.cpp
 *
 *  Created on: Sep 16, 2018
 *      Author: elrond16
 */

#include "array.h"
#include <iostream>
#include "PerformanceTest.h"

using namespace std;

int main() {
    auto t = PerformanceTest<simple_array<double>>(4);
    t.randomize();
    t.printArrays();

    return 0;
}
