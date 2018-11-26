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

struct vector6d {
    double ex = 0, ey = 0, ez = 0, bx = 0, by = 0, bz = 0;
};

template <typename T>
void randomize(T & a) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++){
                a(i, j, k) = gen();
            }
        }
    }
}

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
void randomize_vector6d(T & a) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++){
                a(i, j, k).ex = gen();
                a(i, j, k).ey = gen();
                a(i, j, k).ez = gen();
                a(i, j, k).bx = gen();
                a(i, j, k).by = gen();
                a(i, j, k).bz = gen();
            }
        }
    }
}

template <typename T>
void advance_b(const T & ex, const T & ey, const T & ez, T & bx, T & by, T & bz) {
    const int n = ex.get_n();

    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                bx(i,j,k) += 0.5*(-dtdy*(ez(i,j+1,k)-ez(i,j,k)) + dtdz*(ey(i,j,k+1)-ey(i,j,k)));
                by(i,j,k) += 0.5*( dtdx*(ez(i+1,j,k)-ez(i,j,k)) - dtdz*(ex(i,j,k+1)-ex(i,j,k)));
                bz(i,j,k) += 0.5*( dtdy*(ex(i,j+1,k)-ex(i,j,k)) - dtdx*(ey(i+1,j,k)-ey(i,j,k)));
            }
        }
    }
    {int i=0;
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                by(i,j,k) += 0.5*( dtdx*(ez(i+1,j,k)-ez(i,j,k)) - dtdz*(ex(i,j,k+1)-ex(i,j,k)));
                bz(i,j,k) += 0.5*( dtdy*(ex(i,j+1,k)-ex(i,j,k)) - dtdx*(ey(i+1,j,k)-ey(i,j,k)));
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        {int j=0;
            for(int k=1;k<n-1;k++) {
                bx(i,j,k) += 0.5*(-dtdy*(ez(i,j+1,k)-ez(i,j,k)) + dtdz*(ey(i,j,k+1)-ey(i,j,k)));
                bz(i,j,k) += 0.5*( dtdy*(ex(i,j+1,k)-ex(i,j,k)) - dtdx*(ey(i+1,j,k)-ey(i,j,k)));
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            {int k=0;
                bx(i,j,k) += 0.5*(-dtdy*(ez(i,j+1,k)-ez(i,j,k)) + dtdz*(ey(i,j,k+1)-ey(i,j,k)));
                by(i,j,k) += 0.5*( dtdx*(ez(i+1,j,k)-ez(i,j,k)) - dtdz*(ex(i,j,k+1)-ex(i,j,k)));
            }
        }
    }
}

void advance_b_z(const morton_array<double> & ex, const morton_array<double> & ey, const morton_array<double> & ez,
        morton_array<double> & bx, morton_array<double> & by, morton_array<double> & bz) {
    const unsigned int n = ex.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        if (ex.is_xmax(curr) || ex.is_ymax(curr) || ex.is_zmax(curr)) {
            continue;
        }

        auto x_next = ex.get_x_next(curr);
        auto y_next = ex.get_y_next(curr);
        auto z_next = ex.get_z_next(curr);

        if (ex.is_xmin(curr)) {
            if (ex.is_ymin(curr) || ex.is_zmin(curr)) {
                continue;
            } else {
                by[curr] += 0.5*( dtdx*(ez[x_next]-ez[curr]) - dtdz*(ex[z_next]-ex[curr]));
                bz[curr] += 0.5*( dtdy*(ex[y_next]-ex[curr]) - dtdx*(ey[x_next]-ey[curr]));
            }
        } else if (ex.is_ymin(curr)) {
            if (ex.is_zmin(curr)) {
                continue;
            } else {
                bx[curr] += 0.5*(-dtdy*(ez[y_next]-ez[curr]) + dtdz*(ey[z_next]-ey[curr]));
                bz[curr] += 0.5*( dtdy*(ex[y_next]-ex[curr]) - dtdx*(ey[x_next]-ey[curr]));
            }
        } else if (ex.is_zmin(curr)) {
            bx[curr] += 0.5*(-dtdy*(ez[y_next]-ez[curr]) + dtdz*(ey[z_next]-ey[curr]));
            by[curr] += 0.5*( dtdx*(ez[x_next]-ez[curr]) - dtdz*(ex[z_next]-ex[curr]));
        } else {
            bx[curr] += 0.5*(-dtdy*(ez[y_next]-ez[curr]) + dtdz*(ey[z_next]-ey[curr]));
            by[curr] += 0.5*( dtdx*(ez[x_next]-ez[curr]) - dtdz*(ex[z_next]-ex[curr]));
            bz[curr] += 0.5*( dtdy*(ex[y_next]-ex[curr]) - dtdx*(ey[x_next]-ey[curr]));
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
void advance_b_vector6d(const T & v) {
    const int n = v.get_n();

    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                v(i,j,k).bx += 0.5*(-dtdy*(v(i,j+1,k).ez-v(i,j,k).ez) + dtdz*(v(i,j,k+1).ey-v(i,j,k).ey));
                v(i,j,k).by += 0.5*( dtdx*(v(i+1,j,k).ez-v(i,j,k).ez) - dtdz*(v(i,j,k+1).ex-v(i,j,k).ex));
                v(i,j,k).bz += 0.5*( dtdy*(v(i,j+1,k).ex-v(i,j,k).ex) - dtdx*(v(i+1,j,k).ey-v(i,j,k).ey));
            }
        }
    }
    {int i=0;
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                v(i,j,k).by += 0.5*( dtdx*(v(i+1,j,k).ez-v(i,j,k).ez) - dtdz*(v(i,j,k+1).ex-v(i,j,k).ex));
                v(i,j,k).bz += 0.5*( dtdy*(v(i,j+1,k).ex-v(i,j,k).ex) - dtdx*(v(i+1,j,k).ey-v(i,j,k).ey));
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        {int j=0;
            for(int k=1;k<n-1;k++) {
                v(i,j,k).bx += 0.5*(-dtdy*(v(i,j+1,k).ez-v(i,j,k).ez) + dtdz*(v(i,j,k+1).ey-v(i,j,k).ey));
                v(i,j,k).bz += 0.5*( dtdy*(v(i,j+1,k).ex-v(i,j,k).ex) - dtdx*(v(i+1,j,k).ey-v(i,j,k).ey));
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            {int k=0;
                v(i,j,k).bx += 0.5*(-dtdy*(v(i,j+1,k).ez-v(i,j,k).ez) + dtdz*(v(i,j,k+1).ey-v(i,j,k).ey));
                v(i,j,k).by += 0.5*( dtdx*(v(i+1,j,k).ez-v(i,j,k).ez) - dtdz*(v(i,j,k+1).ex-v(i,j,k).ex));
            }
        }
    }
}

void advance_b_vector6d_z(const morton_array<vector6d> & v) {
    const unsigned int n = v.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        if (v.is_xmax(curr) || v.is_ymax(curr) || v.is_zmax(curr)) {
            continue;
        }

        auto x_next = v.get_x_next(curr);
        auto y_next = v.get_y_next(curr);
        auto z_next = v.get_z_next(curr);

        if (v.is_xmin(curr)) {
            if (v.is_ymin(curr) || v.is_zmin(curr)) {
                continue;
            } else {
                v[curr].by += 0.5*( dtdx*(v[x_next].ez-v[curr].ez) - dtdz*(v[z_next].ex-v[curr].ex));
                v[curr].bz += 0.5*( dtdy*(v[y_next].ex-v[curr].ex) - dtdx*(v[x_next].ey-v[curr].ey));
            }
        } else if (v.is_ymin(curr)) {
            if (v.is_zmin(curr)) {
                continue;
            } else {
                v[curr].bx += 0.5*(-dtdy*(v[y_next].ez-v[curr].ez) + dtdz*(v[z_next].ey-v[curr].ey));
                v[curr].bz += 0.5*( dtdy*(v[y_next].ex-v[curr].ex) - dtdx*(v[x_next].ey-v[curr].ey));
            }
        } else if (v.is_zmin(curr)) {
            v[curr].bx += 0.5*(-dtdy*(v[y_next].ez-v[curr].ez) + dtdz*(v[z_next].ey-v[curr].ey));
            v[curr].by += 0.5*( dtdx*(v[x_next].ez-v[curr].ez) - dtdz*(v[z_next].ex-v[curr].ex));
        } else {
            v[curr].bx += 0.5*(-dtdy*(v[y_next].ez-v[curr].ez) + dtdz*(v[z_next].ey-v[curr].ey));
            v[curr].by += 0.5*( dtdx*(v[x_next].ez-v[curr].ez) - dtdz*(v[z_next].ex-v[curr].ex));
            v[curr].bz += 0.5*( dtdy*(v[y_next].ex-v[curr].ex) - dtdx*(v[x_next].ey-v[curr].ey));
        }
    }
}

template <typename T>
void advance_e(T & ex, T & ey, T & ez, const T & bx, const T & by, const T & bz) {
    const int n = ex.get_n();

    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                ex(i,j,k) +=  dtdy*(bz(i,j,k)-bz(i,j-1,k)) - dtdz*(by(i,j,k)-by(i,j,k-1));
                ey(i,j,k) += -dtdx*(bz(i,j,k)-bz(i-1,j,k)) + dtdz*(bx(i,j,k)-bx(i,j,k-1));
                ez(i,j,k) +=  dtdx*(by(i,j,k)-by(i-1,j,k)) - dtdy*(bx(i,j,k)-bx(i,j-1,k));
            }
        }
    }
    {int i=0;
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                ex(i,j,k) +=  dtdy*(bz(i,j,k)-bz(i,j-1,k)) - dtdz*(by(i,j,k)-by(i,j,k-1));
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        {int j=0;
            for(int k=1;k<n-1;k++) {
                ey(i,j,k) += -dtdx*(bz(i,j,k)-bz(i-1,j,k)) + dtdz*(bx(i,j,k)-bx(i,j,k-1));
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            {int k=0;
                ez(i,j,k) +=  dtdx*(by(i,j,k)-by(i-1,j,k)) - dtdy*(bx(i,j,k)-bx(i,j-1,k));
            }
        }
    }
}

void advance_e_z(morton_array<double> & ex, morton_array<double> & ey, morton_array<double> & ez,
        const morton_array<double> & bx, const morton_array<double> & by, const morton_array<double> & bz) {
    const unsigned int n = ex.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        if (ex.is_xmax(curr) || ex.is_ymax(curr) || ex.is_zmax(curr)) {
            continue;
        }

        auto x_prev = ex.get_x_prev(curr);
        auto y_prev = ex.get_y_prev(curr);
        auto z_prev = ex.get_z_prev(curr);

        if (ex.is_xmin(curr)) {
            if (ex.is_ymin(curr) || ex.is_zmin(curr)) {
                continue;
            } else {
                ex[curr] +=  dtdy*(bz[curr]-bz[y_prev]) - dtdz*(by[curr]-by[z_prev]);
            }
        } else if (ex.is_ymin(curr)) {
            if (ex.is_zmin(curr)) {
                continue;
            } else {
                ey[curr] += -dtdx*(bz[curr]-bz[x_prev]) + dtdz*(bx[curr]-bx[z_prev]);
            }
        } else if (ex.is_zmin(curr)) {
            ez[curr] +=  dtdx*(by[curr]-by[x_prev]) - dtdy*(bx[curr]-bx[y_prev]);
        } else {
            ex[curr] +=  dtdy*(bz[curr]-bz[y_prev]) - dtdz*(by[curr]-by[z_prev]);
            ey[curr] += -dtdx*(bz[curr]-bz[x_prev]) + dtdz*(bx[curr]-bx[z_prev]);
            ez[curr] +=  dtdx*(by[curr]-by[x_prev]) - dtdy*(bx[curr]-bx[y_prev]);
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
void advance_e_vector6d(T & v) {
    const int n = v.get_n();

    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                v(i,j,k).ex +=  dtdy*(v(i,j,k).bz-v(i,j-1,k).bz) - dtdz*(v(i,j,k).by-v(i,j,k-1).by);
                v(i,j,k).ey += -dtdx*(v(i,j,k).bz-v(i-1,j,k).bz) + dtdz*(v(i,j,k).bx-v(i,j,k-1).bx);
                v(i,j,k).ez +=  dtdx*(v(i,j,k).by-v(i-1,j,k).by) - dtdy*(v(i,j,k).bx-v(i,j-1,k).bx);
            }
        }
    }
    {int i=0;
        for(int j=1;j<n-1;j++) {
            for(int k=1;k<n-1;k++) {
                v(i,j,k).ex +=  dtdy*(v(i,j,k).bz-v(i,j-1,k).bz) - dtdz*(v(i,j,k).by-v(i,j,k-1).by);
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        {int j=0;
            for(int k=1;k<n-1;k++) {
                v(i,j,k).ey += -dtdx*(v(i,j,k).bz-v(i-1,j,k).bz) + dtdz*(v(i,j,k).bx-v(i,j,k-1).bx);
            }
        }
    }
    for(int i=1;i<n-1;i++) {
        for(int j=1;j<n-1;j++) {
            int k=0;
            v(i,j,k).ez +=  dtdx*(v(i,j,k).by-v(i-1,j,k).by) - dtdy*(v(i,j,k).bx-v(i,j-1,k).bx);
        }
    }
}

void advance_e_vector6d_z(morton_array<vector6d> & v) {
    const unsigned int n = v.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        if (v.is_xmax(curr) || v.is_ymax(curr) || v.is_zmax(curr)) {
            continue;
        }

        auto x_prev = v.get_x_prev(curr);
        auto y_prev = v.get_y_prev(curr);
        auto z_prev = v.get_z_prev(curr);

        if (v.is_xmin(curr)) {
            if (v.is_ymin(curr) || v.is_zmin(curr)) {
                continue;
            } else {
                v[curr].ex +=  dtdy*(v[curr].bz-v[y_prev].bz) - dtdz*(v[curr].by-v[z_prev].by);
            }
        } else if (v.is_ymin(curr)) {
            if (v.is_zmin(curr)) {
                continue;
            } else {
                v[curr].ey += -dtdx*(v[curr].bz-v[x_prev].bz) + dtdz*(v[curr].bx-v[z_prev].bx);
            }
        } else if (v.is_zmin(curr)) {
            v[curr].ez +=  dtdx*(v[curr].by-v[x_prev].by) - dtdy*(v[curr].bx-v[y_prev].bx);
        } else {
            v[curr].ex +=  dtdy*(v[curr].bz-v[y_prev].bz) - dtdz*(v[curr].by-v[z_prev].by);
            v[curr].ey += -dtdx*(v[curr].bz-v[x_prev].bz) + dtdz*(v[curr].bx-v[z_prev].bx);
            v[curr].ez +=  dtdx*(v[curr].by-v[x_prev].by) - dtdy*(v[curr].bx-v[y_prev].bx);
        }
    }
}

template <typename T>
void advance(T & ex, T & ey, T & ez, T & bx, T & by, T & bz) {
    advance_b(ex, ey, ez, bx, by, bz);
    advance_e(ex, ey, ez, bx, by, bz);
    advance_b(ex, ey, ez, bx, by, bz);
}

void advance_z(morton_array<double> & ex, morton_array<double> & ey, morton_array<double> & ez,
        morton_array<double> & bx, morton_array<double> & by, morton_array<double> & bz) {
    advance_b_z(ex, ey, ez, bx, by, bz);
    advance_e_z(ex, ey, ez, bx, by, bz);
    advance_b_z(ex, ey, ez, bx, by, bz);
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
void advance_vector6d(T & v) {
    advance_b_vector6d(v);
    advance_e_vector6d(v);
    advance_b_vector6d(v);
}

void advance_vector6d_z(morton_array<vector6d> & v) {
    advance_b_vector6d_z(v);
    advance_e_vector6d_z(v);
    advance_b_vector6d_z(v);
}

template <typename T>
void run_test(std::function<void(T&, T&, T&, T&, T&, T&)> func, const std::string & testname, int size, int iterations=10) {
    T ex(size);
    T ey(size);
    T ez(size);
    T bx(size);
    T by(size);
    T bz(size);
    vector<double> times(iterations);
    for (int i = 0; i < iterations; i++) {
        randomize(ex);
        randomize(ey);
        randomize(ez);
        randomize(bx);
        randomize(by);
        randomize(bz);
        MPI_Barrier(MPI_COMM_WORLD);
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        func(ex, ey, ez, bx, by, bz);
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
        double throughput = (6 * static_cast<double>(mpi_size) * size * size * size * sizeof(double)) / GB * 1e3 / avg; // Gb / s
        cout << boost::format("%15s ") % testname << boost::format("%5d: ") % size << boost::format("%9.3f ms") % avg
                << boost::format(" (+-%7.3f ms)") % standard_deviation << boost::format(" %5.2f Gb/s") % throughput
                << endl;
    }
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

template <typename T>
void run_test_vector6d(std::function<void(T&)> func, const std::string & testname, int size, int iterations=10) {
    T a(size);
    vector<double> times(iterations);
    for (int i = 0; i < iterations; i++) {
        randomize_vector6d(a);
        MPI_Barrier(MPI_COMM_WORLD);
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        func(a);
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
        double throughput = (static_cast<double>(mpi_size) * size * size * size * sizeof(vector6d)) / GB * 1e3 / avg; // Gb / s
        cout << boost::format("%15s ") % testname << boost::format("%5d: ") % size << boost::format("%9.3f ms") % avg
                << boost::format(" (+-%7.3f ms)") % standard_deviation << boost::format(" %5.2f Gb/s") % throughput
                << endl;
    }
}

template <typename T, typename V>
bool is_equal(T & a, V & b) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                if (std::fabs(b(i, j, k) - a(i, j, k)) > 1e-6) return false;
            }
        }
    }
    return true;
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

template <typename T, typename V>
bool is_equal_vector6d(T & a, V & b) {
    int n = a.get_n();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                if (std::fabs(b(i, j, k).ex - a(i, j, k).ex) > 1e-6) return false;
                if (std::fabs(b(i, j, k).ey - a(i, j, k).ey) > 1e-6) return false;
                if (std::fabs(b(i, j, k).ez - a(i, j, k).ez) > 1e-6) return false;
                if (std::fabs(b(i, j, k).bx - a(i, j, k).bx) > 1e-6) return false;
                if (std::fabs(b(i, j, k).by - a(i, j, k).by) > 1e-6) return false;
                if (std::fabs(b(i, j, k).bz - a(i, j, k).bz) > 1e-6) return false;
            }
        }
    }
    return true;
}

void run_checks(int size) {
    cout << "Checking size " << size << "..." << endl;

    {
        // ordinary checks

        morton_array<double> ex1(size);
        morton_array<double> ey1(size);
        morton_array<double> ez1(size);
        morton_array<double> bx1(size);
        morton_array<double> by1(size);
        morton_array<double> bz1(size);

        morton_array<double> ex2(size);
        morton_array<double> ey2(size);
        morton_array<double> ez2(size);
        morton_array<double> bx2(size);
        morton_array<double> by2(size);
        morton_array<double> bz2(size);

        randomize(ex1);
        randomize(ey1);
        randomize(ez1);

        ex2 = ex1;
        ey2 = ey1;
        ez2 = ez1;

        randomize(bx1);
        randomize(by1);
        randomize(bz1);

        bx2 = bx1;
        by2 = by1;
        bz2 = bz1;

        advance_b(ex1, ey1, ez1, bx1, by1, bz1);
        advance_b_z(ex2, ey2, ez2, bx2, by2, bz2);

        if (!is_equal(bx1, bx2) || !is_equal(by1, by2) || !is_equal(bz1, bz2))
            cout << "advance_b_z is different" << endl;

        randomize(ex1);
        randomize(ey1);
        randomize(ez1);

        ex2 = ex1;
        ey2 = ey1;
        ez2 = ez1;

        randomize(bx1);
        randomize(by1);
        randomize(bz1);

        bx2 = bx1;
        by2 = by1;
        bz2 = bz1;

        advance_e(ex1, ey1, ez1, bx1, by1, bz1);
        advance_e_z(ex2, ey2, ez2, bx2, by2, bz2);

        if (!is_equal(ex1, ex2) || !is_equal(ey1, ey2) || !is_equal(ez1, ez2))
            cout << "advance_e_z is different" << endl;

        randomize(ex1);
        randomize(ey1);
        randomize(ez1);

        ex2 = ex1;
        ey2 = ey1;
        ez2 = ez1;

        randomize(bx1);
        randomize(by1);
        randomize(bz1);

        bx2 = bx1;
        by2 = by1;
        bz2 = bz1;

        advance(ex1, ey1, ez1, bx1, by1, bz1);
        advance_z(ex2, ey2, ez2, bx2, by2, bz2);

        if (!is_equal(ex1, ex2) || !is_equal(ey1, ey2) || !is_equal(ez1, ez2)
                || !is_equal(bx1, bx2) || !is_equal(by1, by2) || !is_equal(bz1, bz2))
            cout << "advance_z is different" << endl;
    }

    {
        // vector3d checks

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

    {
        // vector6d checks

        morton_array<vector6d> v1(size);
        morton_array<vector6d> v2(size);

        randomize_vector6d(v1);
        v2 = v1;

        advance_b_vector6d(v1);
        advance_b_vector6d_z(v2);

        if (!is_equal_vector6d(v1, v2)) cout << "advance_b_vector6d_z is different" << endl;

        randomize_vector6d(v1);
        v2 = v1;

        advance_e_vector6d(v1);
        advance_e_vector6d_z(v2);

        if (!is_equal_vector6d(v1, v2)) cout << "advance_e_vector6d_z is different" << endl;

        randomize_vector6d(v1);
        v2 = v1;

        advance_vector6d(v1);
        advance_vector6d_z(v2);

        if (!is_equal_vector6d(v1, v2)) cout << "advance_vector6d_z is different" << endl;
    }

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
            cout << endl;
        }

        if (mpi_rank == 0) {
            cout << "scalar" << endl;
        }
        run_test<simple_array<double>>(advance_b<simple_array<double>>, "simple(b)", size, iterations);
        run_test<cached_array<double>>(advance_b<cached_array<double>>, "cached(b)", size, iterations);
        run_test<morton_array<double>>(advance_b<morton_array<double>>, "morton_bad(b)", size, iterations);
        run_test<morton_array<double>>(advance_b_z, "morton(b)", size, iterations);

        run_test<simple_array<double>>(advance_e<simple_array<double>>, "simple(e)", size, iterations);
        run_test<cached_array<double>>(advance_e<cached_array<double>>, "cached(e)", size, iterations);
        run_test<morton_array<double>>(advance_e<morton_array<double>>, "morton_bad(e)", size, iterations);
        run_test<morton_array<double>>(advance_e_z, "morton(e)", size, iterations);

        run_test<simple_array<double>>(advance<simple_array<double>>, "simple", size, iterations);
        run_test<cached_array<double>>(advance<cached_array<double>>, "cached", size, iterations);
        run_test<morton_array<double>>(advance<morton_array<double>>, "morton_bad", size, iterations);
        run_test<morton_array<double>>(advance_z, "morton", size, iterations);

        if (mpi_rank == 0) {
            cout << "vector 3d" << endl;
        }
        run_test_vector<simple_array<vector3d>>(advance_b_vector<simple_array<vector3d>>, "simple(b)", size, iterations);
        run_test_vector<cached_array<vector3d>>(advance_b_vector<cached_array<vector3d>>, "cached(b)", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_b_vector<morton_array<vector3d>>, "morton_bad(b)", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_b_vector_z, "morton(b)", size, iterations);

        run_test_vector<simple_array<vector3d>>(advance_e_vector<simple_array<vector3d>>, "simple(e)", size, iterations);
        run_test_vector<cached_array<vector3d>>(advance_e_vector<cached_array<vector3d>>, "cached(e)", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_e_vector<morton_array<vector3d>>, "morton_bad(e)", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_e_vector_z, "morton(e)", size, iterations);

        run_test_vector<simple_array<vector3d>>(advance_vector<simple_array<vector3d>>, "simple", size, iterations);
        run_test_vector<cached_array<vector3d>>(advance_vector<cached_array<vector3d>>, "cached", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_vector<morton_array<vector3d>>, "morton_bad", size, iterations);
        run_test_vector<morton_array<vector3d>>(advance_vector_z, "morton", size, iterations);

        if (mpi_rank == 0) {
            cout << "vector 6d" << endl;
        }

        run_test_vector6d<simple_array<vector6d>>(advance_b_vector6d<simple_array<vector6d>>, "simple(b)", size, iterations);
        run_test_vector6d<cached_array<vector6d>>(advance_b_vector6d<cached_array<vector6d>>, "cached(b)", size, iterations);
        run_test_vector6d<morton_array<vector6d>>(advance_b_vector6d<morton_array<vector6d>>, "morton_bad(b)", size, iterations);
        run_test_vector6d<morton_array<vector6d>>(advance_b_vector6d_z, "morton(b)", size, iterations);

        run_test_vector6d<simple_array<vector6d>>(advance_e_vector6d<simple_array<vector6d>>, "simple(e)", size, iterations);
        run_test_vector6d<cached_array<vector6d>>(advance_e_vector6d<cached_array<vector6d>>, "cached(e)", size, iterations);
        run_test_vector6d<morton_array<vector6d>>(advance_e_vector6d<morton_array<vector6d>>, "morton_bad(e)", size, iterations);
        run_test_vector6d<morton_array<vector6d>>(advance_e_vector6d_z, "morton(e)", size, iterations);

        run_test_vector6d<simple_array<vector6d>>(advance_vector6d<simple_array<vector6d>>, "simple", size, iterations);
        run_test_vector6d<cached_array<vector6d>>(advance_vector6d<cached_array<vector6d>>, "cached", size, iterations);
        run_test_vector6d<morton_array<vector6d>>(advance_vector6d<morton_array<vector6d>>, "morton_bad", size, iterations);
        run_test_vector6d<morton_array<vector6d>>(advance_vector6d_z, "morton", size, iterations);
    }

    MPI_Finalize();
    return 0;
}

