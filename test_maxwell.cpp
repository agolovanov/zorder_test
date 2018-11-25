#include <random>
#include <functional>
#include <vector>
#include <chrono>
#include <boost/format.hpp>
#include <algorithm>
#include "mpi.h"
#include "array3d.h"

using namespace std;

uniform_real_distribution<double> generator { -10, 10 };
default_random_engine re;

int mpi_rank;
int mpi_size;

const double GB = 1024 * 1024 * 1024;

double dt;
double dtdx;
double dtdy;
double dtdz;

double gen() {
    return generator(re);
}

struct vector3d {
    double x = 0, y = 0, z = 0;
};

template <typename T>
void randomize_vector(T & a) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++){
                a(i, j, k).x = gen();
                a(i, j, k).y = gen();
                a(i, j, k).z = gen();
            }
        }
    }
}

template <typename T>
void advance_b_vector(const T & e, T & b) {
    const int n = e.get_n();

    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                b(i,j,k).x += 0.5*(-dtdy*(e(i,j+1,k).z-e(i,j,k).z) + dtdz*(e(i,j,k+1).y-e(i,j,k).y));
                b(i,j,k).y += 0.5*( dtdx*(e(i+1,j,k).z-e(i,j,k).z) - dtdz*(e(i,j,k+1).x-e(i,j,k).x));
                b(i,j,k).z += 0.5*( dtdy*(e(i,j+1,k).x-e(i,j,k).x) - dtdx*(e(i+1,j,k).y-e(i,j,k).y));
            }
        }
    }
    {int i=0;
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                b(i,j,k).y += 0.5*( dtdx*(e(i+1,j,k).z-e(i,j,k).z) - dtdz*(e(i,j,k+1).x-e(i,j,k).x));
                b(i,j,k).z += 0.5*( dtdy*(e(i,j+1,k).x-e(i,j,k).x) - dtdx*(e(i+1,j,k).y-e(i,j,k).y));
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        {int j=0;
            for(int k=1;k<n-1;k++) {
                b(i,j,k).x += 0.5*(-dtdy*(e(i,j+1,k).z-e(i,j,k).z) + dtdz*(e(i,j,k+1).y-e(i,j,k).y));
                b(i,j,k).z += 0.5*( dtdy*(e(i,j+1,k).x-e(i,j,k).x) - dtdx*(e(i+1,j,k).y-e(i,j,k).y));
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            {int k=0;
                b(i,j,k).x += 0.5*(-dtdy*(e(i,j+1,k).z-e(i,j,k).z) + dtdz*(e(i,j,k+1).y-e(i,j,k).y));
                b(i,j,k).y += 0.5*( dtdx*(e(i+1,j,k).z-e(i,j,k).z) - dtdz*(e(i,j,k+1).x-e(i,j,k).x));
            }
        }
    }
}

void advance_b_vector_z(const morton_array<vector3d> & e, morton_array<vector3d> & b) {
    const unsigned int n = e.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        if (e.is_xmax(curr) || e.is_ymax(curr) || e.is_zmax(curr)) {
            continue;
        }

        auto x_next = e.get_x_next(curr);
        auto y_next = e.get_y_next(curr);
        auto z_next = e.get_z_next(curr);

        if (e.is_xmin(curr)) {
            if (e.is_ymin(curr) || e.is_zmin(curr)) {
                continue;
            } else {
                b[curr].y += 0.5*( dtdx*(e[x_next].z-e[curr].z) - dtdz*(e[z_next].x-e[curr].x));
                b[curr].z += 0.5*( dtdy*(e[y_next].x-e[curr].x) - dtdx*(e[x_next].y-e[curr].y));
            }
        } else if (e.is_ymin(curr)) {
            if (e.is_zmin(curr)) {
                continue;
            } else {
                b[curr].x += 0.5*(-dtdy*(e[y_next].z-e[curr].z) + dtdz*(e[z_next].y-e[curr].y));
                b[curr].z += 0.5*( dtdy*(e[y_next].x-e[curr].x) - dtdx*(e[x_next].y-e[curr].y));
            }
        } else if (e.is_zmin(curr)) {
            b[curr].x += 0.5*(-dtdy*(e[y_next].z-e[curr].z) + dtdz*(e[z_next].y-e[curr].y));
            b[curr].y += 0.5*( dtdx*(e[x_next].z-e[curr].z) - dtdz*(e[z_next].x-e[curr].x));
        } else {
            b[curr].x += 0.5*(-dtdy*(e[y_next].z-e[curr].z) + dtdz*(e[z_next].y-e[curr].y));
            b[curr].y += 0.5*( dtdx*(e[x_next].z-e[curr].z) - dtdz*(e[z_next].x-e[curr].x));
            b[curr].z += 0.5*( dtdy*(e[y_next].x-e[curr].x) - dtdx*(e[x_next].y-e[curr].y));
        }
    }
}

template <typename T>
void advance_e_vector(T & e, const T & b) {
    const int n = e.get_n();

    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                e(i,j,k).x +=  dtdy*(b(i,j,k).z-b(i,j-1,k).z) - dtdz*(b(i,j,k).y-b(i,j,k-1).y);
                e(i,j,k).y += -dtdx*(b(i,j,k).z-b(i-1,j,k).z) + dtdz*(b(i,j,k).x-b(i,j,k-1).x);
                e(i,j,k).z +=  dtdx*(b(i,j,k).y-b(i-1,j,k).y) - dtdy*(b(i,j,k).x-b(i,j-1,k).x);
            }
        }
    }
    {int i=0;
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                e(i,j,k).x +=  dtdy*(b(i,j,k).z-b(i,j-1,k).z) - dtdz*(b(i,j,k).y-b(i,j,k-1).y);
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        {int j=0;
            for(int k=1;k<n-1;k++) {
                e(i,j,k).y += -dtdx*(b(i,j,k).z-b(i-1,j,k).z) + dtdz*(b(i,j,k).x-b(i,j,k-1).x);
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            {int k=0;
                e(i,j,k).z +=  dtdx*(b(i,j,k).y-b(i-1,j,k).y) - dtdy*(b(i,j,k).x-b(i,j-1,k).x);
            }
        }
    }
}

void advance_e_vector_z(morton_array<vector3d> & e, const morton_array<vector3d> & b) {
    const unsigned int n = e.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        if (e.is_xmax(curr) || e.is_ymax(curr) || e.is_zmax(curr)) {
            continue;
        }

        auto x_prev = e.get_x_prev(curr);
        auto y_prev = e.get_y_prev(curr);
        auto z_prev = e.get_z_prev(curr);

        if (e.is_xmin(curr)) {
            if (e.is_ymin(curr) || e.is_zmin(curr)) {
                continue;
            } else {
                e[curr].x +=  dtdy*(b[curr].z-b[y_prev].z) - dtdz*(b[curr].y-b[z_prev].y);
            }
        } else if (e.is_ymin(curr)) {
            if (e.is_zmin(curr)) {
                continue;
            } else {
                e[curr].y += -dtdx*(b[curr].z-b[x_prev].z) + dtdz*(b[curr].x-b[z_prev].x);
            }
        } else if (e.is_zmin(curr)) {
            e[curr].z +=  dtdx*(b[curr].y-b[x_prev].y) - dtdy*(b[curr].x-b[y_prev].x);
        } else {
            e[curr].x +=  dtdy*(b[curr].z-b[y_prev].z) - dtdz*(b[curr].y-b[z_prev].y);
            e[curr].y += -dtdx*(b[curr].z-b[x_prev].z) + dtdz*(b[curr].x-b[z_prev].x);
            e[curr].z +=  dtdx*(b[curr].y-b[x_prev].y) - dtdy*(b[curr].x-b[y_prev].x);
        }
    }
}

template <typename T>
void advance_vector(T & e, T & b) {
    advance_b_vector(e, b);
    advance_e_vector(e, b);
    advance_b_vector(e, b);
}

void advance_vector_z(morton_array<vector3d> & e, morton_array<vector3d> & b) {
    advance_b_vector_z(e, b);
    advance_e_vector_z(e, b);
    advance_b_vector_z(e, b);
}

template <typename T>
void run_test_vector(std::function<void(T&, T&)> func, const std::string & testname, int size, int iterations=10) {
    T a(size);
    T b(size);
    vector<double> times(iterations);
    for (int i = 0; i < iterations; i++) {
        randomize_vector(a);
        randomize_vector(b);
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
        double throughput = (2 * static_cast<double>(mpi_size) * size * size * size * sizeof(vector3d)) / GB * 1e3 / avg; // Gb / s
        cout << boost::format("%15s ") % testname << boost::format("%5d: ") % size << boost::format("%9.3f ms") % avg
                << boost::format(" (+-%7.3f ms)") % standard_deviation << boost::format(" %5.2f Gb/s") % throughput
                << endl;
    }
}

template <typename T, typename V>
bool is_equal_vector(T & a, V & b) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                if (std::fabs(b(i, j, k).x - a(i, j, k).x) > 1e-6) return false;
                if (std::fabs(b(i, j, k).y - a(i, j, k).y) > 1e-6) return false;
                if (std::fabs(b(i, j, k).z - a(i, j, k).z) > 1e-6) return false;
            }
        }
    }
    return true;
}

void run_checks(int size) {
    cout << "Checking size " << size << "..." << endl;

    morton_array<vector3d> e1(size);
    morton_array<vector3d> e2(size);
    morton_array<vector3d> b1(size);
    morton_array<vector3d> b2(size);

    randomize_vector(e1);
    e2 = e1;
    randomize_vector(b1);
    b2 = b1;

    advance_b_vector(e1, b1);
    advance_b_vector_z(e2, b2);

    if (!is_equal_vector(b1, b2)) cout << "advance_b_vector_z is different" << endl;

    randomize_vector(e1);
    e2 = e1;
    randomize_vector(b1);
    b2 = b1;

    advance_e_vector(e1, b1);
    advance_e_vector_z(e2, b2);

    if (!is_equal_vector(e1, e2)) cout << "advance_e_vector_z is different" << endl;

    randomize_vector(e1);
    e2 = e1;
    randomize_vector(b1);
    b2 = b1;

    advance_vector(e1, b1);
    advance_vector_z(e2, b2);

    if (!is_equal_vector(e1, e2) || !is_equal_vector(b1, b2)) cout << "advance_vector_z is different" << endl;
}

template <typename T>
void print_array_vector(T& a) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                cout << boost::format("[%6.2f,%6.2f,%6.2f]") % a(i, j, k).x % a(i, j, k).y % a(i, j, k).z;
            }
            cout << endl;
        }
        cout << endl;
    }
}

int main(int argc, char **argv) {
    dt = 0.01;
    double dx = 0.05;
    double dy = 0.05;
    double dz = 0.05;

    dtdx = dt / dx;
    dtdy = dt / dy;
    dtdz = dt / dz;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    re.seed(chrono::system_clock::now().time_since_epoch().count());

    const int iterations = 6;

    vector<int> sizes = {32, 64, 128, 256, 512};

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

    for (auto size : sizes) {
        if (mpi_rank == 0) {
            cout << "Size " << size << endl;
        }

        run_test_vector<simple_array<vector3d>>(advance_b_vector<simple_array<vector3d>>, "s_vector(b)", size, iterations);
        run_test_vector<cached_array<vector3d>>(advance_b_vector<cached_array<vector3d>>, "c_vector(b)", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_b_vector<morton_array<vector3d>>, "m_vector(b)", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_b_vector_z, "m_vector_z(b)", size, iterations);

        run_test_vector<simple_array<vector3d>>(advance_e_vector<simple_array<vector3d>>, "s_vector(e)", size, iterations);
        run_test_vector<cached_array<vector3d>>(advance_e_vector<cached_array<vector3d>>, "c_vector(e)", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_e_vector<morton_array<vector3d>>, "m_vector(e)", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_e_vector_z, "m_vector_z(e)", size, iterations);

        run_test_vector<simple_array<vector3d>>(advance_vector<simple_array<vector3d>>, "s_vector", size, iterations);
        run_test_vector<cached_array<vector3d>>(advance_vector<cached_array<vector3d>>, "c_vector", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_vector<morton_array<vector3d>>, "m_vector", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_vector_z, "m_vector_z", size, iterations);
    }

    MPI_Finalize();
    return 0;
}

