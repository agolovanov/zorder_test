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
#include <omp.h>
#include "mpi.h"
#include <algorithm>

using namespace std;

uniform_real_distribution<double> generator { -10, 10 };
default_random_engine re;

const double GB = 1024 * 1024 * 1024;

int mpi_rank;
int mpi_size;

double gen() {
    return generator(re);
}

template <typename T>
void randomize(T & a) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a(i, j) = gen();
        }
    }
}

template<typename T>
void print_array(T& a) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << boost::format("%6.2f") % a(i, j);
        }
        cout << endl;
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
        MPI_Barrier(MPI_COMM_WORLD);
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        func(a, b);
        MPI_Barrier(MPI_COMM_WORLD);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        times[i] = chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1e6;
    }
    if (mpi_rank == 0) {
        sort(times.begin(), times.end());
        times.pop_back();
        double avg = accumulate(times.begin(), times.end(), 0.0) / times.size();
        double avg_sq = inner_product(times.begin(), times.end(), times.begin(), 0.0) / times.size();
        double standard_deviation = sqrt(avg_sq - avg * avg);
        double throughput = (static_cast<double>(mpi_size) * size * size * sizeof(double)) / GB * 1e3 / avg; // Gb / s
        cout << boost::format("%15s ") % testname << boost::format("%5d: ") % size << boost::format("%9.3f ms") % avg
                << boost::format(" (+-%7.3f ms)") % standard_deviation << boost::format(" %5.2f Gb/s") % throughput
                << endl;
    }
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
    int n = a.get_n();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b(i, j) = a(i, j);
        }
    }
}

template <typename T>
void copy_wo_test(T & a, T & b) {
    int n = a.get_n();
    #pragma omp parallel for
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            b(i, j) = a(i, j);
        }
    }
}

void copy_z(morton_array<double> & a, morton_array<double> & b) {
    int n = a.get_size();
    for (int i = 0; i < n; i++) {
        b[i] = a[i];
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
        b(0, j) = a(0, j) + a(1, j);
    }
    #pragma omp parallel for
    for (int i = 1; i < n - 1; i++) {
        for (int j = 0; j < n; j++) {
            b(i, j) = a(i-1, j) + a(i, j) + a(i+1, j);
        }
    }
    for (int j = 0; j < n; j++) {
        b(n-1, j) = a(n-2, j) + a(n-1, j);
    }
}

void sum_triplets_bad_z(morton_array<double> & a, morton_array<double> & b) {
    int n = a.get_size();
    for (int i = 0; i < n; i++) {
        if (a.is_ymin(i)) {
            b[i] = a[i] + a[a.get_y_next(i)];
        } else if (a.is_ymax(i)) {
            b[i] = a[i] + a[a.get_y_prev(i)];
        } else {
            b[i] = a[i] + a[a.get_y_next(i)] + a[a.get_y_prev(i)];
        }
    }
}

template <typename T>
void sum_triplets_good(T & a, T & b) {
    int n = a.get_n();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        b(i, 0) = a(i, 0) + a(i, 1);
        for (int j = 1; j < n-1; j++) {
            b(i, j) = a(i, j-1) + a(i, j) + a(i, j+1);
        }
        b(i, n-1) = a(i, n-1) + a(i, n-2);
    }
}

void sum_triplets_good_z(morton_array<double> & a, morton_array<double> & b) {
    int n = a.get_size();
    for (int i = 0; i < n; i++) {
        if (a.is_xmin(i)) {
            b[i] = a[i] + a[a.get_x_next(i)];
        } else if (a.is_xmax(i)) {
            b[i] = a[i] + a[a.get_x_prev(i)];
        } else {
            b[i] = a[i] + a[a.get_x_next(i)] + a[a.get_x_prev(i)];
        }
    }
}

template <typename T>
void sum_pentlets_bad(T & a, T & b) {
    int n = a.get_n();
    for (int j = 0; j < n; j++) {
        b(0, j) = a(0, j) + a(1, j) + a(2, j);
    }
    for (int j = 0; j < n; j++) {
        b(1, j) = a(0, j) + a(1, j) + a(2, j) + a(3, j);
    }
    #pragma omp parallel for
    for (int i = 2; i < n - 2; i++) {
        for (int j = 0; j < n; j++) {
            b(i, j) = a(i-2, j) + a(i-1, j) + a(i, j) + a(i+1, j) + a(i+2, j);
        }
    }
    for (int j = 0; j < n; j++) {
        b(n-2, j) = a(n-4, j) + a(n-3, j) + a(n-2, j) + a(n-1, j);
    }
    for (int j = 0; j < n; j++) {
        b(n-1, j) = a(n-3, j) + a(n-2, j) + a(n-1, j);
    }
}

void sum_pentlets_bad_z(morton_array<double> & a, morton_array<double> & b) {
    int n = a.get_size();
    for (int i = 0; i < n; i++) {
        if (a.is_ymin(i)) {
            int i1 = a.get_y_next(i);
            int i2 = a.get_y_next(i1);
            b[i] = a[i] + a[i1] + a[i2];
        } else if (a.is_ymax(i)) {
            int i1 = a.get_y_prev(i);
            int i2 = a.get_y_prev(i1);
            b[i] = a[i] + a[i1] + a[i2];
        } else {
            int i1_prev = a.get_y_prev(i);
            int i1_next = a.get_y_next(i);
            if (a.is_ymin(i1_prev)) {
                int i2_next = a.get_y_next(i1_next);
                b[i] = a[i1_prev] + a[i] + a[i1_next] + a[i2_next];
            } else if (a.is_ymax(i1_next)) {
                int i2_prev = a.get_y_prev(i1_prev);
                b[i] = a[i2_prev] + a[i1_prev] + a[i] + a[i1_next];
            } else {
                int i2_next = a.get_y_next(i1_next);
                int i2_prev = a.get_y_prev(i1_prev);
                b[i] = a[i2_prev] + a[i1_prev] + a[i] + a[i1_next] + a[i2_next];
            }
        }
    }
}

template <typename T>
void sum_pentlets_good(T & a, T & b) {
    int n = a.get_n();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        b(i, 0) = a(i, 0) + a(i, 1) + a(i, 2);
        b(i, 1) = a(i, 0) + a(i, 1) + a(i, 2) + a(i, 3);
        for (int j = 2; j < n-2; j++) {
            b(i, j) = a(i, j-2) + a(i, j-1) + a(i, j) + a(i, j+1) + a(i, j+2);
        }
        b(i, n-2) = a(i, n-4) + a(i, n-3) + a(i, n-2) + a(i, n-1);
        b(i, n-1) = a(i, n-3) + a(i, n-2) + a(i, n-1);
    }
}

void sum_pentlets_good_z(morton_array<double> & a, morton_array<double> & b) {
    int n = a.get_size();
    for (int i = 0; i < n; i++) {
        if (a.is_xmin(i)) {
            int i1 = a.get_x_next(i);
            int i2 = a.get_x_next(i1);
            b[i] = a[i] + a[i1] + a[i2];
        } else if (a.is_xmax(i)) {
            int i1 = a.get_x_prev(i);
            int i2 = a.get_x_prev(i1);
            b[i] = a[i] + a[i1] + a[i2];
        } else {
            int i1_prev = a.get_x_prev(i);
            int i1_next = a.get_x_next(i);
            if (a.is_xmin(i1_prev)) {
                int i2_next = a.get_x_next(i1_next);
                b[i] = a[i1_prev] + a[i] + a[i1_next] + a[i2_next];
            } else if (a.is_xmax(i1_next)) {
                int i2_prev = a.get_x_prev(i1_prev);
                b[i] = a[i2_prev] + a[i1_prev] + a[i] + a[i1_next];
            } else {
                int i2_next = a.get_x_next(i1_next);
                int i2_prev = a.get_x_prev(i1_prev);
                b[i] = a[i2_prev] + a[i1_prev] + a[i] + a[i1_next] + a[i2_next];
            }
        }
    }
}

template <typename T>
void sum_neighbors(T & a, T & b) {
    int n = a.get_n();
    #pragma omp parallel for
    for (int i = 0; i < n-1; i++) {
        for (int j = 0; j < n-1; j++) {
            b(i, j) = a(i, j) + a(i+1, j) + a(i, j+1) + a(i+1, j+1);
        }
        b(i, n-1) = a(i, n-1) + a(i+1, n-1);
    }
    for (int j = 0; j < n-1; j++) {
        b(n-1, j) = a(n-1, j) + a(n-1, j+1);
    }
    b(n-1, n-1) = a(n-1, n-1);
}

void sum_neighbors_z(morton_array<double> & a, morton_array<double> & b) {
    int n = a.get_size();
    for (int i = 0; i < n; i++) {
        if (a.is_xmax(i)) {
            if (a.is_ymax(i)) {
                b[i] = a[i];
            } else {
                b[i] = a[i] + a[a.get_y_next(i)];
            }
        } else if (a.is_ymax(i)) {
            b[i] = a[i] + a[a.get_x_next(i)];
        } else {
            int x_next = a.get_x_next(i);
            int y_next = a.get_y_next(i);
            int xy_next = a.get_x_next(y_next);
            b[i] = a[i] + a[x_next] + a[y_next] + a[xy_next];
        }
    }
}

template <typename T, typename V>
void copy(T & a, V & b) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b(i, j) = a(i, j);
        }
    }
}

template <typename T, typename V>
bool is_equal(T & a, V & b) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (std::fabs(b(i, j) - a(i, j)) > 1e-6) return false;
        }
    }
    return true;
}

void run_checks(int size) {
    cout << "Checking size " << size << "..." << endl;
    morton_array<double> a(size);

    {
        cached_array<double> b(size);
        simple_array<double> c(size);

        randomize(a);
        copy(a, b);
        copy(a, c);
        if (!is_equal(a, b)) cout << "a and b are different" << endl;
        if (!is_equal(a, c)) cout << "a and c are different" << endl;
        if (!is_equal(b, c)) cout << "b and c are different" << endl;
    }

    morton_array<double> b(size);
    morton_array<double> c(size);

    copy(a, b);
    copy_test(a, c);
    if (!is_equal(b, c)) cout << "copy_test is different" << endl;

    copy(a, b);
    copy_wo_test(a, c);
    if (!is_equal(b, c)) cout << "copy_wo_test is different" << endl;

    copy(a, b);
    copy_true_test(a, c);
    if (!is_equal(b, c)) cout << "copy_true_test is different" << endl;

    copy(a, b);
    copy_z(a, c);
    if (!is_equal(b, c)) cout << "copy_z is different" << endl;

    sum_triplets_bad(a, b);
    sum_triplets_bad_z(a, c);
    if (!is_equal(b, c)) cout << "sum_triplets_bad_z is different" << endl;

    sum_triplets_good(a, b);
    sum_triplets_good_z(a, c);
    if (!is_equal(b, c)) cout << "sum_triplets_good_z is different" << endl;

    sum_pentlets_bad(a, b);
    sum_pentlets_bad_z(a, c);
    if (!is_equal(b, c)) cout << "sum_pentlets_bad_z is different" << endl;

    sum_pentlets_good(a, b);
    sum_pentlets_good_z(a, c);
    if (!is_equal(b, c)) cout << "sum_pentlets_good_z is different" << endl;

    sum_neighbors(a, b);
    sum_neighbors_z(a, c);
    if (!is_equal(b, c)) cout << "sum_neighbors_z is different" << endl;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    re.seed(chrono::system_clock::now().time_since_epoch().count());

    const int iterations = 6;

    vector<int> sizes = {512, 1024, 2048, 4096, 8192};

    if (argc > 1) {
        sizes.erase(sizes.begin(), sizes.end());
        for (int i = 1; i < argc; i++) {
            sizes.push_back(atoi(argv[i]));
        }
    }

    if (mpi_rank == 0) {
        cout << "Sizes:";
        for (auto size : sizes) {
            cout << " " << size;
        }
        cout << endl;
    }

    if (mpi_rank == 0) {
        for (auto size : sizes) {
            run_checks(size);
        }
    }

    if (mpi_rank == 0) {
        cout << "Morton array:" << endl;
    }

    run_test<morton_array<double>>(copy_test<morton_array<double>>, "copy", sizes, iterations);
    run_test<morton_array<double>>(copy_wo_test<morton_array<double>>, "copy_wo", sizes, iterations);
    run_test<morton_array<double>>(copy_z, "copy_z", sizes, iterations);
    run_test<morton_array<double>>(copy_true_test<morton_array<double>>, "copy_true", sizes, iterations);
    run_test<morton_array<double>>(sum_triplets_bad<morton_array<double>>, "triplets_bad", sizes, iterations);
    run_test<morton_array<double>>(sum_triplets_bad_z, "triplets_bad_z", sizes, iterations);
    run_test<morton_array<double>>(sum_triplets_good<morton_array<double>>, "triplets_good", sizes, iterations);
    run_test<morton_array<double>>(sum_triplets_good_z, "triplets_good_z", sizes, iterations);
    run_test<morton_array<double>>(sum_pentlets_bad<morton_array<double>>, "pentlets_bad", sizes, iterations);
    run_test<morton_array<double>>(sum_pentlets_bad_z, "pentlets_bad_z", sizes, iterations);
    run_test<morton_array<double>>(sum_pentlets_good<morton_array<double>>, "pentlets_good", sizes, iterations);
    run_test<morton_array<double>>(sum_pentlets_good_z, "pentlets_good_z", sizes, iterations);
    run_test<morton_array<double>>(sum_neighbors<morton_array<double>>, "sum_neighbors", sizes, iterations);
    run_test<morton_array<double>>(sum_neighbors_z, "sum_neighbors_z", sizes, iterations);

    if (mpi_rank == 0) {
        cout << "Simple array:" << endl;
    }

    run_test<simple_array<double>>(copy_test<simple_array<double>>, "copy", sizes, iterations);
    run_test<simple_array<double>>(copy_wo_test<simple_array<double>>, "copy_wo", sizes, iterations);
    run_test<simple_array<double>>(copy_true_test<simple_array<double>>, "copy_true", sizes, iterations);
    run_test<simple_array<double>>(sum_triplets_bad<simple_array<double>>, "triplets_bad", sizes, iterations);
    run_test<simple_array<double>>(sum_triplets_good<simple_array<double>>, "triplets_good", sizes, iterations);
    run_test<simple_array<double>>(sum_pentlets_bad<simple_array<double>>, "pentlets_bad", sizes, iterations);
    run_test<simple_array<double>>(sum_pentlets_good<simple_array<double>>, "pentlets_good", sizes, iterations);
    run_test<simple_array<double>>(sum_neighbors<simple_array<double>>, "sum_neighbors", sizes, iterations);

    if (mpi_rank == 0) {
        cout << "Cached array:" << endl;
    }

    run_test<cached_array<double>>(copy_test<cached_array<double>>, "copy", sizes, iterations);
    run_test<cached_array<double>>(copy_wo_test<cached_array<double>>, "copy_wo", sizes, iterations);
    run_test<cached_array<double>>(copy_true_test<cached_array<double>>, "copy_true", sizes, iterations);
    run_test<cached_array<double>>(sum_triplets_bad<cached_array<double>>, "triplets_bad", sizes, iterations);
    run_test<cached_array<double>>(sum_triplets_good<cached_array<double>>, "triplets_good", sizes, iterations);
    run_test<cached_array<double>>(sum_pentlets_bad<cached_array<double>>, "pentlets_bad", sizes, iterations);
    run_test<cached_array<double>>(sum_pentlets_good<cached_array<double>>, "pentlets_good", sizes, iterations);
    run_test<cached_array<double>>(sum_neighbors<cached_array<double>>, "sum_neighbors", sizes, iterations);

    MPI_Finalize();
    return 0;
}
