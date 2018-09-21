/*
 * main.cpp
 *
 *  Created on: Sep 16, 2018
 *      Author: elrond16
 */

#include "array.h"
#include <iostream>
#include <random>
#include <functional>
#include <chrono>
#include <boost/format.hpp>
#include <vector>
#include <cmath>

using namespace std;

uniform_real_distribution<double> generator { -10, 10 };
default_random_engine re;

const double GB = 1024 * 1024 * 1024;

double gen() {
    return generator(re);
}

template <typename T>
void randomize(T & a) {
    int n1 = a.get_n1();
    int n2 = a.get_n2();
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            a(i, j) = gen();
        }
    }
}

template<typename T>
void print_array(T& a) {
    int n1 = a.get_n1();
    int n2 = a.get_n2();
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            cout << boost::format("%6.2f") % a(i, j);
        }
        cout << endl;
    }
}

template <typename T>
void run_test(std::function<void(T&, T&)> func, const std::string & testname, int size, int iterations=10) {
    T a(size, size);
    T b(size, size);
    randomize(b);
    vector<double> times(iterations);
    for (int i = 0; i < iterations; i++) {
        randomize(a);
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        func(a, b);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        times[i] = chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1e6;
    }
    double avg = accumulate(times.begin(), times.end(), 0.0) / times.size();
    double avg_sq = inner_product(times.begin(), times.end(), times.begin(), 0.0) / times.size();
    double standard_deviation = sqrt(avg_sq - avg * avg);
    double throughput = (size * size * sizeof(double)) / GB * 1e3 / avg; // Gb / s
    cout << boost::format("%15s ") % testname << boost::format("%5d: ") % size << boost::format("%9.3f ms") % avg
            << boost::format(" (+-%7.3f ms)") % standard_deviation << boost::format(" %4.2f Gb/s") % throughput << endl;
}

template <typename T>
void run_test(std::function<void(T&, T&)> func, const std::string & testname, const vector<int> & sizes, int iterations=10) {
    for (int size : sizes) {
        run_test(func, testname, size, iterations);
    }
}

template <typename T>
void empty_test(T & a, T & b) { }

template <typename T>
void copy_test(T & a, T & b) {
    int n1 = a.get_n1();
    int n2 = a.get_n2();
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            b(i, j) = a(i, j);
        }
    }
}

template <typename T>
void copy_wo_test(T & a, T & b) {
    int n1 = a.get_n1();
    int n2 = a.get_n2();
    for (int j = 0; j < n2; j++) {
        for (int i = 0; i < n1; i++) {
            b(i, j) = a(i, j);
        }
    }
}

template <typename T>
void copy_true_test(T & a, T & b) {
    b = a;
}

template <typename T>
void sum_triplets_bad(T & a, T & b) {
    int n1 = a.get_n1();
    int n2 = a.get_n2();
    for (int j = 0; j < n2; j++) {
        b(0, j) = a(0, j) + a(1, j);
    }
    for (int i = 1; i < n1 - 1; i++) {
        for (int j = 0; j < n2; j++) {
            b(i, j) = a(i-1, j) + a(i, j) + a(i+1, j);
        }
    }
    for (int j = 0; j < n2; j++) {
        b(n1-1, j) = a(n1-2, j) + a(n1-1, j);
    }
}

template <typename T>
void sum_triplets_good(T & a, T & b) {
    int n1 = a.get_n1();
    int n2 = a.get_n2();
    for (int i = 0; i < n1; i++) {
        b(i, 0) = a(i, 0) + a(i, 1);
        for (int j = 1; j < n2-1; j++) {
            b(i, j) = a(i, j-1) + a(i, j) + a(i, j+1);
        }
        b(i, n2-1) = a(i, n2-1) + a(i, n2-2);
    }
}

template <typename T>
void sum_pentlets_bad(T & a, T & b) {
    int n1 = a.get_n1();
    int n2 = a.get_n2();
    for (int j = 0; j < n2; j++) {
        b(0, j) = a(0, j) + a(1, j) + a(2, j);
    }
    for (int j = 0; j < n2; j++) {
        b(1, j) = a(0, j) + a(1, j) + a(2, j) + a(3, j);
    }
    for (int i = 2; i < n1 - 2; i++) {
        for (int j = 0; j < n2; j++) {
            b(i, j) = a(i-2, j) + a(i-1, j) + a(i, j) + a(i+1, j) + a(i+2, j);
        }
    }
    for (int j = 0; j < n2; j++) {
        b(n1-2, j) = a(n1-4, j) + a(n1-3, j) + a(n1-2, j) + a(n1-1, j);
    }
    for (int j = 0; j < n2; j++) {
        b(n1-1, j) = a(n1-3, j) + a(n1-2, j) + a(n1-1, j);
    }
}

template <typename T>
void sum_pentlets_good(T & a, T & b) {
    int n1 = a.get_n1();
    int n2 = a.get_n2();
    for (int i = 0; i < n1; i++) {
        b(i, 0) = a(i, 0) + a(i, 1) + a(i, 2);
        b(i, 1) = a(i, 0) + a(i, 1) + a(i, 2) + a(i, 3);
        for (int j = 2; j < n2-2; j++) {
            b(i, j) = a(i, j-2) + a(i, j-1) + a(i, j) + a(i, j+1) + a(i, j+2);
        }
        b(i, n2-2) = a(i, n2-4) + a(i, n2-3) + a(i, n2-2) + a(i, n2-1);
        b(i, n2-1) = a(i, n2-3) + a(i, n2-2) + a(i, n2-1);
    }
}

int main() {
    /*
    simple_array<double> a(4, 4);
    randomize(a);
    print_array(a);
    simple_array<double> b(4, 4);
    copy_true_test(a, b);
    print_array(b);

    return 0;
    */
    const int iterations = 5;

    vector<int> sizes = {512, 1024, 2048, 4096, 8192};

    cout << "Simple array:" << endl;

    //run_test<simple_array<double>>(empty_test<simple_array<double>>, "empty", sizes, 10);
    run_test<simple_array<double>>(copy_test<simple_array<double>>, "copy", sizes, iterations);
    run_test<simple_array<double>>(copy_wo_test<simple_array<double>>, "copy_wo", sizes, iterations);
    run_test<simple_array<double>>(copy_true_test<simple_array<double>>, "copy_true", sizes, iterations);
    run_test<simple_array<double>>(sum_triplets_bad<simple_array<double>>, "triplets_bad", sizes, iterations);
    run_test<simple_array<double>>(sum_triplets_good<simple_array<double>>, "triplets_good", sizes, iterations);
    run_test<simple_array<double>>(sum_pentlets_bad<simple_array<double>>, "pentlets_bad", sizes, iterations);
    run_test<simple_array<double>>(sum_pentlets_good<simple_array<double>>, "pentlets_good", sizes, iterations);

    cout << "Cached array:" << endl;

    run_test<cached_array<double>>(copy_test<cached_array<double>>, "copy", sizes, iterations);
    run_test<cached_array<double>>(copy_wo_test<cached_array<double>>, "copy_wo", sizes, iterations);
    run_test<cached_array<double>>(copy_true_test<cached_array<double>>, "copy_true", sizes, iterations);
    run_test<cached_array<double>>(sum_triplets_bad<cached_array<double>>, "triplets_bad", sizes, iterations);
    run_test<cached_array<double>>(sum_triplets_good<cached_array<double>>, "triplets_good", sizes, iterations);
    run_test<cached_array<double>>(sum_pentlets_bad<cached_array<double>>, "pentlets_bad", sizes, iterations);
    run_test<cached_array<double>>(sum_pentlets_good<cached_array<double>>, "pentlets_good", sizes, iterations);

    return 0;
}
