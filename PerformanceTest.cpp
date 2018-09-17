/*
 * PerformanceTest.cpp
 *
 *  Created on: Sep 17, 2018
 *      Author: agolovanov
 */

#include "PerformanceTest.h"
#include "array.h"
#include <random>
#include <iostream>
#include <boost/format.hpp>
#include <vector>
#include <chrono>

template<typename T>
PerformanceTest<T>::PerformanceTest(int size, int iterations) :
        size(size), iterations(iterations), a(size, size), b(size, size), c(size, size) {
}

template<typename T>
void PerformanceTest<T>::randomize() {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            a(i, j) = gen();
            b(i, j) = gen();
            c(i, j) = gen();
        }
    }
}

template<typename T>
void PerformanceTest<T>::printArray(T& a) {
    using namespace std;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << boost::format("%6.2f") % a(i, j);
        }
        cout << endl;
    }
}

template<typename T>
void PerformanceTest<T>::printArrays() {
    using namespace std;
    cout << "a:" << endl;
    printArray(a);
    cout << "b:" << endl;
    printArray(b);
    cout << "c:" << endl;
    printArray(c);
}

template<typename T>
void PerformanceTest<T>::runTest(std::function<void(T&, T&, T&)> func, std::string testname) {
    std::vector<double> times(iterations);
    for (int i = 0; i < iterations; i++) {
        randomize();
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        func(a, b, c);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        times[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    }
}

template class PerformanceTest<cached_array<double>> ;
template class PerformanceTest<simple_array<double>> ;
