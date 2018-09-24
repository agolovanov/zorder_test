#include "array3d.h"
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
    int n = a.get_n();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++){
                a(i, j, k) = gen();
            }
        }
    }
}

template <typename T>
void run_test(std::function<void(T&, T&)> func, const std::string & testname, int size, int iterations=10) {
    T a(size);
    T b(size);
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
    double throughput = (size * size * size * sizeof(double)) / GB * 1e3 / avg; // Gb / s
    cout << boost::format("%15s ") % testname << boost::format("%5d: ") % size << boost::format("%9.3f ms") % avg
            << boost::format(" (+-%7.3f ms)") % standard_deviation << boost::format(" %5.2f Gb/s") % throughput << endl;
}

template <typename T>
void run_test(std::function<void(T&, T&)> func, const std::string & testname, const vector<int> & sizes, int iterations=10) {
    for (int size : sizes) {
        run_test(func, testname, size, iterations);
    }
}

template <typename T>
void copy_test(T & a, T & b) {
    int n = a.get_n();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                b(i, j, k) = a(i, j, k);
            }
        }
    }
}

template <typename T>
void copy_wo_test(T & a, T & b) {
    int n = a.get_n();
    #pragma omp parallel for
    for (int k = 0; k < n; k++) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                b(i, j, k) = a(i, j, k);
            }
        }
    }
}

template <typename T>
void copy_true_test(T & a, T & b) {
    b = a;
}

template <typename T>
void sum_triplets_bad(T & a, T & b) {
    int n = a.get_n();
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            b(0, j, k) = a(0, j, k) + a(1, j, k);
        }
    }
    #pragma omp parallel for
    for (int i = 1; i < n - 1; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                b(i, j, k) = a(i-1, j, k) + a(i, j, k) + a(i+1, j, k);
            }
        }
    }
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            b(n-1, j, k) = a(n-2, j, k) + a(n-1, j, k);
        }
    }
}

template <typename T>
void sum_triplets_good(T & a, T & b) {
    int n = a.get_n();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b(i, j, 0) = a(i, j, 0) + a(i, j, 1);
            for (int k = 1; k < n-1; k++) {
                b(i, j, k) = a(i, j, k-1) + a(i, j, k) + a(i, j, k+1);
            }
            b(i, j, n-1) = a(i, j, n-1) + a(i, j, n-2);
        }
    }
}


template <typename T>
void sum_pentlets_bad(T & a, T & b) {
    int n = a.get_n();
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            b(0, j, k) = a(0, j, k) + a(1, j, k) + a(2, j, k);
        }
    }
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            b(1, j, k) = a(0, j, k) + a(1, j, k) + a(2, j, k) + a(3, j, k);
        }
    }
    #pragma omp parallel for
    for (int i = 2; i < n - 2; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                b(i, j, k) = a(i-2, j, k) + a(i-1, j, k) + a(i, j, k) + a(i+1, j, k) + a(i+2, j, k);
            }
        }
    }
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            b(n-2, j, k) = a(n-4, j, k) + a(n-3, j, k) + a(n-2, j, k) + a(n-1, j, k);
        }
    }
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            b(n-1, j, k) = a(n-3, j, k) + a(n-2, j, k) + a(n-1, j, k);
        }
    }
}

template <typename T>
void sum_pentlets_good(T & a, T & b) {
    int n = a.get_n();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b(i, j, 0) = a(i, j, 0) + a(i, j, 1) + a(i, j, 2);
            b(i, j, 1) = a(i, j, 0) + a(i, j, 1) + a(i, j, 2) + a(i, j, 3);
            for (int k = 2; k < n-2; k++) {
                b(i, j, k) = a(i, j, k-2) + a(i, j, k-1) + a(i, j, k) + a(i, j, k+1) + a(i, j, k+2);
            }
            b(i, j, n-2) = a(i, j, n-4) + a(i, j, n-3) + a(i, j, n-2) + a(i, j, n-1);
            b(i, j, n-1) = a(i, j, n-3) + a(i, j, n-2) + a(i, j, n-1);
        }
    }
}


int main() {
    const int iterations = 5;

    vector<int> sizes = {32, 64, 128, 256, 512};

    cout << "Morton array:" << endl;

    run_test<morton_array<double>>(copy_test<morton_array<double>>, "copy", sizes, iterations);
    run_test<morton_array<double>>(copy_wo_test<morton_array<double>>, "copy_wo", sizes, iterations);
    run_test<morton_array<double>>(copy_true_test<morton_array<double>>, "copy_true", sizes, iterations);
    run_test<morton_array<double>>(sum_triplets_bad<morton_array<double>>, "triplets_bad", sizes, iterations);
    run_test<morton_array<double>>(sum_triplets_good<morton_array<double>>, "triplets_good", sizes, iterations);
    run_test<morton_array<double>>(sum_pentlets_bad<morton_array<double>>, "pentlets_bad", sizes, iterations);
    run_test<morton_array<double>>(sum_pentlets_good<morton_array<double>>, "pentlets_good", sizes, iterations);

    cout << "Simple array:" << endl;

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
