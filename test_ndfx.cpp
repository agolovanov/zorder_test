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

const double kappa = 1.25 * (sqrt(3) - 1) * 0.5 / sqrt(3);
const double bx = (1 - kappa);
const double ax = 0.5 * (1 - bx);

double dt;
double dtdx;
double dtdy;
double dtdz;
double ay, az;
double by, bz;

double gen() {
    return generator(re);
}

struct vector3d {
    double x = 0, y = 0, z = 0;

    vector3d operator-(const vector3d & other) {
        vector3d res;
        res.x = x - other.x;
        res.y = y - other.y;
        res.z = z - other.z;
        return res;
    }
};

template <typename T>
void randomize(T & a) {
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
void advance_b(const T & e, T & b) {
    const int n = e.get_n();

    for(int i=2;i<n-1;i++) {
        for(int j=2;j<n-1;j++) {
            for(int k=2;k<n-1;k++) {
                b(i,j,k).x += 0.5*(   -dtdy*( bz*(e(i,j+1,k).z-e(i,j,k).z) + az*(e(i,j+1,k+1).z-e(i,j,k+1).z+e(i,j+1,k-1).z-e(i,j,k-1).z) )    +   dtdz*( by*(e(i,j,k+1).y-e(i,j,k).y) + ay*(e(i,j+1,k+1).y-e(i,j+1,k).y+e(i,j-1,k+1).y-e(i,j-1,k).y) )  );
                b(i,j,k).y += 0.5*(   dtdx*( bz*(e(i+1,j,k).z-e(i,j,k).z) + az*(e(i+1,j,k+1).z-e(i,j,k+1).z+e(i+1,j,k-1).z-e(i,j,k-1).z) )    -   dtdz*( bx*(e(i,j,k+1).x-e(i,j,k).x) + ax*(e(i+1,j,k+1).x-e(i+1,j,k).x+e(i-1,j,k+1).x-e(i-1,j,k).x) )  );
                b(i,j,k).z += 0.5*(   dtdy*( bx*(e(i,j+1,k).x-e(i,j,k).x) + ax*(e(i+1,j+1,k).x-e(i+1,j,k).x+e(i-1,j+1,k).x-e(i-1,j,k).x) )    -   dtdx*( by*(e(i+1,j,k).y-e(i,j,k).y) + ay*(e(i+1,j+1,k).y-e(i,j+1,k).y+e(i+1,j-1,k).y-e(i,j-1,k).y) )  );
            }
        }
    }
    {int i=1;
        for(int j=2;j<n-1;j++) {
            for(int k=2;k<n-1;k++) {
                b(i,j,k).y += 0.5*(   dtdx*( bz*(e(i+1,j,k).z-e(i,j,k).z) + az*(e(i+1,j,k+1).z-e(i,j,k+1).z+e(i+1,j,k-1).z-e(i,j,k-1).z) )    -   dtdz*( bx*(e(i,j,k+1).x-e(i,j,k).x) + ax*(e(i+1,j,k+1).x-e(i+1,j,k).x+e(i-1,j,k+1).x-e(i-1,j,k).x) )  );
                b(i,j,k).z += 0.5*(   dtdy*( bx*(e(i,j+1,k).x-e(i,j,k).x) + ax*(e(i+1,j+1,k).x-e(i+1,j,k).x+e(i-1,j+1,k).x-e(i-1,j,k).x) )    -   dtdx*( by*(e(i+1,j,k).y-e(i,j,k).y) + ay*(e(i+1,j+1,k).y-e(i,j+1,k).y+e(i+1,j-1,k).y-e(i,j-1,k).y) )  );
            }
        }
    }
    for(int i=2;i<n-1;i++) {
        {int j=1;
            for(int k=2;k<n-1;k++) {
                b(i,j,k).x += 0.5*(   -dtdy*( bz*(e(i,j+1,k).z-e(i,j,k).z) + az*(e(i,j+1,k+1).z-e(i,j,k+1).z+e(i,j+1,k-1).z-e(i,j,k-1).z) )    +   dtdz*( by*(e(i,j,k+1).y-e(i,j,k).y) + ay*(e(i,j+1,k+1).y-e(i,j+1,k).y+e(i,j-1,k+1).y-e(i,j-1,k).y) )  );
                b(i,j,k).z += 0.5*(   dtdy*( bx*(e(i,j+1,k).x-e(i,j,k).x) + ax*(e(i+1,j+1,k).x-e(i+1,j,k).x+e(i-1,j+1,k).x-e(i-1,j,k).x) )    -   dtdx*( by*(e(i+1,j,k).y-e(i,j,k).y) + ay*(e(i+1,j+1,k).y-e(i,j+1,k).y+e(i+1,j-1,k).y-e(i,j-1,k).y) )  );
            }
        }
    }
    for(int i=2;i<n-1;i++) {
        for(int j=2;j<n-1;j++) {
            {int k=1;
                b(i,j,k).x += 0.5*(   -dtdy*( bz*(e(i,j+1,k).z-e(i,j,k).z) + az*(e(i,j+1,k+1).z-e(i,j,k+1).z+e(i,j+1,k-1).z-e(i,j,k-1).z) )    +   dtdz*( by*(e(i,j,k+1).y-e(i,j,k).y) + ay*(e(i,j+1,k+1).y-e(i,j+1,k).y+e(i,j-1,k+1).y-e(i,j-1,k).y) )  );
                b(i,j,k).y += 0.5*(   dtdx*( bz*(e(i+1,j,k).z-e(i,j,k).z) + az*(e(i+1,j,k+1).z-e(i,j,k+1).z+e(i+1,j,k-1).z-e(i,j,k-1).z) )    -   dtdz*( bx*(e(i,j,k+1).x-e(i,j,k).x) + ax*(e(i+1,j,k+1).x-e(i+1,j,k).x+e(i-1,j,k+1).x-e(i-1,j,k).x) )  );
            }
        }
    }
}

void advance_b_z(const morton_array<vector3d> & e, morton_array<vector3d> & b) {
    const unsigned int n = e.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        if (e.is_xmax(curr) || e.is_ymax(curr) || e.is_zmax(curr)) {
            continue;
        }

        if (e.is_xmin(curr) || e.is_ymin(curr) || e.is_zmin(curr)) {
            continue;
        }

        auto x_m = e.get_x_prev(curr);
        auto y_m = e.get_y_prev(curr);
        auto z_m = e.get_z_prev(curr);

        auto x_p = e.get_x_next(curr);
        auto y_p = e.get_y_next(curr);
        auto z_p = e.get_z_next(curr);

        auto xy_pp = e.get_y_next(x_p);
        auto xy_pm = e.get_y_prev(x_p);
        auto xy_mp = e.get_y_next(x_m);

        auto xz_pp = e.get_z_next(x_p);
        auto xz_pm = e.get_z_prev(x_p);
        auto xz_mp = e.get_z_next(x_m);

        auto yz_pp = e.get_z_next(y_p);
        auto yz_pm = e.get_z_prev(y_p);
        auto yz_mp = e.get_z_next(y_m);

        if (e.is_xmin(x_m)) {
            if (!e.is_ymin(y_m) && !e.is_zmin(z_m)) {
                b[curr].y += 0.5*(    dtdx*( bz*(e[x_p].z-e[curr].z) + az*(e[xz_pp].z-e[z_p].z+e[xz_pm].z-e[z_m].z) )    -   dtdz*( bx*(e[z_p].x-e[curr].x) + ax*(e[xz_pp].x-e[x_p].x+e[xz_mp].x-e[x_m].x) )  );
                b[curr].z += 0.5*(    dtdy*( bx*(e[y_p].x-e[curr].x) + ax*(e[xy_pp].x-e[x_p].x+e[xy_mp].x-e[x_m].x) )    -   dtdx*( by*(e[x_p].y-e[curr].y) + ay*(e[xy_pp].y-e[y_p].y+e[xy_pm].y-e[y_m].y) )  );
            }
        } else if (e.is_ymin(y_m)) {
            if (!e.is_xmin(x_m) && !e.is_zmin(z_m)) {
                b[curr].x += 0.5*(   -dtdy*( bz*(e[y_p].z-e[curr].z) + az*(e[yz_pp].z-e[z_p].z+e[yz_pm].z-e[z_m].z) )    +   dtdz*( by*(e[z_p].y-e[curr].y) + ay*(e[yz_pp].y-e[y_p].y+e[yz_mp].y-e[y_m].y) )  );
                b[curr].z += 0.5*(    dtdy*( bx*(e[y_p].x-e[curr].x) + ax*(e[xy_pp].x-e[x_p].x+e[xy_mp].x-e[x_m].x) )    -   dtdx*( by*(e[x_p].y-e[curr].y) + ay*(e[xy_pp].y-e[y_p].y+e[xy_pm].y-e[y_m].y) )  );
            }
        } else if (e.is_zmin(z_m)) {
            if (!e.is_xmin(x_m) && !e.is_ymin(y_m)) {
                b[curr].x += 0.5*(   -dtdy*( bz*(e[y_p].z-e[curr].z) + az*(e[yz_pp].z-e[z_p].z+e[yz_pm].z-e[z_m].z) )    +   dtdz*( by*(e[z_p].y-e[curr].y) + ay*(e[yz_pp].y-e[y_p].y+e[yz_mp].y-e[y_m].y) )  );
                b[curr].y += 0.5*(    dtdx*( bz*(e[x_p].z-e[curr].z) + az*(e[xz_pp].z-e[z_p].z+e[xz_pm].z-e[z_m].z) )    -   dtdz*( bx*(e[z_p].x-e[curr].x) + ax*(e[xz_pp].x-e[x_p].x+e[xz_mp].x-e[x_m].x) )  );
            }
        } else {
            b[curr].x += 0.5*(   -dtdy*( bz*(e[y_p].z-e[curr].z) + az*(e[yz_pp].z-e[z_p].z+e[yz_pm].z-e[z_m].z) )    +   dtdz*( by*(e[z_p].y-e[curr].y) + ay*(e[yz_pp].y-e[y_p].y+e[yz_mp].y-e[y_m].y) )  );
            b[curr].y += 0.5*(    dtdx*( bz*(e[x_p].z-e[curr].z) + az*(e[xz_pp].z-e[z_p].z+e[xz_pm].z-e[z_m].z) )    -   dtdz*( bx*(e[z_p].x-e[curr].x) + ax*(e[xz_pp].x-e[x_p].x+e[xz_mp].x-e[x_m].x) )  );
            b[curr].z += 0.5*(    dtdy*( bx*(e[y_p].x-e[curr].x) + ax*(e[xy_pp].x-e[x_p].x+e[xy_mp].x-e[x_m].x) )    -   dtdx*( by*(e[x_p].y-e[curr].y) + ay*(e[xy_pp].y-e[y_p].y+e[xy_pm].y-e[y_m].y) )  );
        }
    }
}

template <typename T>
inline void advance_b_point(const T & e, T & b, int i, int j, int k, int i_prev, int j_prev, int k_prev, int i_next, int j_next, int k_next) {
    b(i,j,k).x += 0.5*(   -dtdy*( bz*(e(i,j_next,k).z-e(i,j,k).z) + az*(e(i,j_next,k_next).z-e(i,j,k_next).z+e(i,j_next,k_prev).z-e(i,j,k_prev).z) )    +   dtdz*( by*(e(i,j,k_next).y-e(i,j,k).y) + ay*(e(i,j_next,k_next).y-e(i,j_next,k).y+e(i,j_prev,k_next).y-e(i,j_prev,k).y) )  );
    b(i,j,k).y += 0.5*(   dtdx*( bz*(e(i_next,j,k).z-e(i,j,k).z) + az*(e(i_next,j,k_next).z-e(i,j,k_next).z+e(i_next,j,k_prev).z-e(i,j,k_prev).z) )    -   dtdz*( bx*(e(i,j,k_next).x-e(i,j,k).x) + ax*(e(i_next,j,k_next).x-e(i_next,j,k).x+e(i_prev,j,k_next).x-e(i_prev,j,k).x) )  );
    b(i,j,k).z += 0.5*(   dtdy*( bx*(e(i,j_next,k).x-e(i,j,k).x) + ax*(e(i_next,j_next,k).x-e(i_next,j,k).x+e(i_prev,j_next,k).x-e(i_prev,j,k).x) )    -   dtdx*( by*(e(i_next,j,k).y-e(i,j,k).y) + ay*(e(i_next,j_next,k).y-e(i,j_next,k).y+e(i_next,j_prev,k).y-e(i,j_prev,k).y) )  );
}

template <typename T>
void advance_b_periodic(const T & e, T & b) {
    const int n = e.get_n();

    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
            for(int k=0;k<n;k++) {
                int i_prev = (i != 0 ? i-1 : n-1);
                int j_prev = (j != 0 ? j-1 : n-1);
                int k_prev = (k != 0 ? k-1 : n-1);
                int i_next = (i != n-1 ? i+1 : 0);
                int j_next = (j != n-1 ? j+1 : 0);
                int k_next = (k != n-1 ? k+1 : 0);
                advance_b_point(e, b, i, j, k, i_prev, j_prev, k_prev, i_next, j_next, k_next);
            }
        }
    }
}

template <typename T>
void advance_b_periodic2(const T & e, T & b) {
    const int n = e.get_n();

    int i = 0;
    int j = 0;
    advance_b_point(e, b, i, j, 0, n-1, n-1, n-1, i+1, j+1, 1);

    for(int k=1;k<n-1;k++) {
        advance_b_point(e, b, i, j, k, n-1, n-1, k-1, i+1, j+1, k+1);
    }

    advance_b_point(e, b, i, j, n-1, n-1, n-1, n-2, i+1, j+1, 0);

    for(int j=1;j<n-1;j++) {
        advance_b_point(e, b, i, j, 0, n-1, j-1, n-1, i+1, j+1, 1);

        for(int k=1;k<n-1;k++) {
            advance_b_point(e, b, i, j, k, n-1, j-1, k-1, i+1, j+1, k+1);
        }

        advance_b_point(e, b, i, j, n-1, n-1, j-1, n-2, i+1, j+1, 0);
    }

    j = n-1;
    advance_b_point(e, b, i, j, 0, n-1, j-1, n-1, i+1, 0, 1);

    for(int k=1;k<n-1;k++) {
        advance_b_point(e, b, i, j, k, n-1, j-1, k-1, i+1, 0, k+1);
    }

    advance_b_point(e, b, i, j, n-1, n-1, j-1, n-2, i+1, 0, 0);

    for(int i=1;i<n-1;i++) {
        int j = 0;
        advance_b_point(e, b, i, j, 0, i-1, n-1, n-1, i+1, j+1, 1);

        for(int k=1;k<n-1;k++) {
            advance_b_point(e, b, i, j, k, i-1, n-1, k-1, i+1, j+1, k+1);
        }

        advance_b_point(e, b, i, j, n-1, i-1, n-1, n-2, i+1, j+1, 0);

        for(int j=1;j<n-1;j++) {
            advance_b_point(e, b, i, j, 0, i-1, j-1, n-1, i+1, j+1, 1);

            for(int k=1;k<n-1;k++) {
                advance_b_point(e, b, i, j, k, i-1, j-1, k-1, i+1, j+1, k+1);
            }

            advance_b_point(e, b, i, j, n-1, i-1, j-1, n-2, i+1, j+1, 0);
        }

        j = n-1;
        advance_b_point(e, b, i, j, 0, i-1, j-1, n-1, i+1, 0, 1);

        for(int k=1;k<n-1;k++) {
            advance_b_point(e, b, i, j, k, i-1, j-1, k-1, i+1, 0, k+1);
        }

        advance_b_point(e, b, i, j, n-1, i-1, j-1, n-2, i+1, 0, 0);
    }

    i = n-1;
    j = 0;
    advance_b_point(e, b, i, j, 0, i-1, n-1, n-1, 0, j+1, 1);

    for(int k=1;k<n-1;k++) {
        advance_b_point(e, b, i, j, k, i-1, n-1, k-1, 0, j+1, k+1);
    }

    advance_b_point(e, b, i, j, n-1, i-1, n-1, n-2, 0, j+1, 0);

    for(int j=1;j<n-1;j++) {
        advance_b_point(e, b, i, j, 0, i-1, j-1, n-1, 0, j+1, 1);

        for(int k=1;k<n-1;k++) {
            advance_b_point(e, b, i, j, k, i-1, j-1, k-1, 0, j+1, k+1);
        }

        advance_b_point(e, b, i, j, n-1, i-1, j-1, n-2, 0, j+1, 0);
    }

    j = n-1;
    advance_b_point(e, b, i, j, 0, i-1, j-1, n-1, 0, 0, 1);

    for(int k=1;k<n-1;k++) {
        advance_b_point(e, b, i, j, k, i-1, j-1, k-1, 0, 0, k+1);
    }

    advance_b_point(e, b, i, j, n-1, i-1, j-1, n-2, 0, 0, 0);
}


void advance_b_z_periodic(const morton_array<vector3d> & e, morton_array<vector3d> & b) {
    const unsigned int n = e.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        auto x_p = e.get_x_next(curr);
        auto y_p = e.get_y_next(curr);
        auto z_p = e.get_z_next(curr);

        auto x_m = e.get_x_prev(curr);
        auto y_m = e.get_y_prev(curr);
        auto z_m = e.get_z_prev(curr);

        auto xy_pp = e.get_y_next(x_p);
        auto xy_pm = e.get_y_prev(x_p);
        auto xy_mp = e.get_y_next(x_m);

        auto xz_pp = e.get_z_next(x_p);
        auto xz_pm = e.get_z_prev(x_p);
        auto xz_mp = e.get_z_next(x_m);

        auto yz_pp = e.get_z_next(y_p);
        auto yz_pm = e.get_z_prev(y_p);
        auto yz_mp = e.get_z_next(y_m);

        b[curr].x += 0.5*(   -dtdy*( bz*(e[y_p].z-e[curr].z) + az*(e[yz_pp].z-e[z_p].z+e[yz_pm].z-e[z_m].z) )    +   dtdz*( by*(e[z_p].y-e[curr].y) + ay*(e[yz_pp].y-e[y_p].y+e[yz_mp].y-e[y_m].y) )  );
        b[curr].y += 0.5*(    dtdx*( bz*(e[x_p].z-e[curr].z) + az*(e[xz_pp].z-e[z_p].z+e[xz_pm].z-e[z_m].z) )    -   dtdz*( bx*(e[z_p].x-e[curr].x) + ax*(e[xz_pp].x-e[x_p].x+e[xz_mp].x-e[x_m].x) )  );
        b[curr].z += 0.5*(    dtdy*( bx*(e[y_p].x-e[curr].x) + ax*(e[xy_pp].x-e[x_p].x+e[xy_mp].x-e[x_m].x) )    -   dtdx*( by*(e[x_p].y-e[curr].y) + ay*(e[xy_pp].y-e[y_p].y+e[xy_pm].y-e[y_m].y) )  );

    }
}

template <typename T>
void advance_e(T & e, const T & b) {
    const int n = e.get_n();

    for(int i=2;i<n-1;i++) {
        for(int j=2;j<n-1;j++) {
            for(int k=2;k<n-1;k++) {
                e(i,j,k).x = e(i,j,k).x   +   dtdy*( bz*(b(i,j,k).z-b(i,j-1,k).z) + az*(b(i,j,k+1).z-b(i,j-1,k+1).z+b(i,j,k-1).z-b(i,j-1,k-1).z) )  -   dtdz*( by*(b(i,j,k).y-b(i,j,k-1).y) + ay*(b(i,j+1,k).y-b(i,j+1,k-1).y+b(i,j-1,k).y-b(i,j-1,k-1).y) );
                e(i,j,k).y = e(i,j,k).y   -   dtdx*( bz*(b(i,j,k).z-b(i-1,j,k).z) + az*(b(i,j,k+1).z-b(i-1,j,k+1).z+b(i,j,k-1).z-b(i-1,j,k-1).z) )  +   dtdz*( bx*(b(i,j,k).x-b(i,j,k-1).x) + ax*(b(i+1,j,k).x-b(i+1,j,k-1).x+b(i-1,j,k).x-b(i-1,j,k-1).x) );
                e(i,j,k).z = e(i,j,k).z   +   dtdx*( by*(b(i,j,k).y-b(i-1,j,k).y) + ay*(b(i,j+1,k).y-b(i-1,j+1,k).y+b(i,j-1,k).y-b(i-1,j-1,k).y) )  -   dtdy*( bx*(b(i,j,k).x-b(i,j-1,k).x) + ax*(b(i+1,j,k).x-b(i+1,j-1,k).x+b(i-1,j,k).x-b(i-1,j-1,k).x) );
            }
        }
    }
    {int i=1;
        for(int j=2;j<n-1;j++) {
            for(int k=2;k<n-1;k++) {
                e(i,j,k).x = e(i,j,k).x   +   dtdy*( bz*(b(i,j,k).z-b(i,j-1,k).z) + az*(b(i,j,k+1).z-b(i,j-1,k+1).z+b(i,j,k-1).z-b(i,j-1,k-1).z) )  -   dtdz*( by*(b(i,j,k).y-b(i,j,k-1).y) + ay*(b(i,j+1,k).y-b(i,j+1,k-1).y+b(i,j-1,k).y-b(i,j-1,k-1).y) );
            }
        }
    }
    for(int i=2;i<n-1;i++) {
        {int j=1;
            for(int k=2;k<n-1;k++) {
                e(i,j,k).y = e(i,j,k).y   -   dtdx*( bz*(b(i,j,k).z-b(i-1,j,k).z) + az*(b(i,j,k+1).z-b(i-1,j,k+1).z+b(i,j,k-1).z-b(i-1,j,k-1).z) )  +   dtdz*( bx*(b(i,j,k).x-b(i,j,k-1).x) + ax*(b(i+1,j,k).x-b(i+1,j,k-1).x+b(i-1,j,k).x-b(i-1,j,k-1).x) );
            }
        }
    }
    for(int i=2;i<n-1;i++) {
        for(int j=2;j<n-1;j++) {
            {int k=1;
                e(i,j,k).z = e(i,j,k).z   +   dtdx*( by*(b(i,j,k).y-b(i-1,j,k).y) + ay*(b(i,j+1,k).y-b(i-1,j+1,k).y+b(i,j-1,k).y-b(i-1,j-1,k).y) )  -   dtdy*( bx*(b(i,j,k).x-b(i,j-1,k).x) + ax*(b(i+1,j,k).x-b(i+1,j-1,k).x+b(i-1,j,k).x-b(i-1,j-1,k).x) );
            }
        }
    }
}

void advance_e_z(morton_array<vector3d> & e, const morton_array<vector3d> & b) {
    const unsigned int n = e.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        if (b.is_xmax(curr) || b.is_ymax(curr) || b.is_zmax(curr)) {
            continue;
        }

        if (b.is_xmin(curr) || b.is_ymin(curr) || b.is_zmin(curr)) {
            continue;
        }

        auto x_p = b.get_x_next(curr);
        auto y_p = b.get_y_next(curr);
        auto z_p = b.get_z_next(curr);

        auto x_m = b.get_x_prev(curr);
        auto y_m = b.get_y_prev(curr);
        auto z_m = b.get_z_prev(curr);

        auto xy_mm = b.get_y_prev(x_m);
        auto xy_pm = b.get_y_prev(x_p);
        auto xy_mp = b.get_y_next(x_m);

        auto xz_mm = b.get_z_prev(x_m);
        auto xz_pm = b.get_z_prev(x_p);
        auto xz_mp = b.get_z_next(x_m);

        auto yz_mm = b.get_z_prev(y_m);
        auto yz_pm = b.get_z_prev(y_p);
        auto yz_mp = b.get_z_next(y_m);

        if (b.is_xmin(x_m)) {
            if (!b.is_ymin(y_m) && !b.is_zmin(z_m)) {
                e[curr].x = e[curr].x   +   dtdy*( bz*(b[curr].z-b[y_m].z) + az*(b[z_p].z-b[yz_mp].z+b[z_m].z-b[yz_mm].z) )  -   dtdz*( by*(b[curr].y-b[z_m].y) + ay*(b[y_p].y-b[yz_pm].y+b[y_m].y-b[yz_mm].y) );
            }
        } else if (b.is_ymin(y_m)) {
            if (!b.is_xmin(x_m) && !b.is_zmin(z_m)) {
                e[curr].y = e[curr].y   -   dtdx*( bz*(b[curr].z-b[x_m].z) + az*(b[z_p].z-b[xz_mp].z+b[z_m].z-b[xz_mm].z) )  +   dtdz*( bx*(b[curr].x-b[z_m].x) + ax*(b[x_p].x-b[xz_pm].x+b[x_m].x-b[xz_mm].x) );            }
        } else if (b.is_zmin(z_m)) {
            if (!b.is_xmin(x_m) && !b.is_ymin(y_m)) {
                e[curr].z = e[curr].z   +   dtdx*( by*(b[curr].y-b[x_m].y) + ay*(b[y_p].y-b[xy_mp].y+b[y_m].y-b[xy_mm].y) )  -   dtdy*( bx*(b[curr].x-b[y_m].x) + ax*(b[x_p].x-b[xy_pm].x+b[x_m].x-b[xy_mm].x) );
            }
        } else {
            e[curr].x = e[curr].x   +   dtdy*( bz*(b[curr].z-b[y_m].z) + az*(b[z_p].z-b[yz_mp].z+b[z_m].z-b[yz_mm].z) )  -   dtdz*( by*(b[curr].y-b[z_m].y) + ay*(b[y_p].y-b[yz_pm].y+b[y_m].y-b[yz_mm].y) );
            e[curr].y = e[curr].y   -   dtdx*( bz*(b[curr].z-b[x_m].z) + az*(b[z_p].z-b[xz_mp].z+b[z_m].z-b[xz_mm].z) )  +   dtdz*( bx*(b[curr].x-b[z_m].x) + ax*(b[x_p].x-b[xz_pm].x+b[x_m].x-b[xz_mm].x) );
            e[curr].z = e[curr].z   +   dtdx*( by*(b[curr].y-b[x_m].y) + ay*(b[y_p].y-b[xy_mp].y+b[y_m].y-b[xy_mm].y) )  -   dtdy*( bx*(b[curr].x-b[y_m].x) + ax*(b[x_p].x-b[xy_pm].x+b[x_m].x-b[xy_mm].x) );
        }

    }
}

template <typename T>
inline void advance_e_point(T & e, const T & b, int i, int j, int k, int i_prev, int j_prev, int k_prev, int i_next, int j_next, int k_next) {
    e(i,j,k).x = e(i,j,k).x   +   dtdy*( bz*(b(i,j,k).z-b(i,j_prev,k).z) + az*(b(i,j,k_next).z-b(i,j_prev,k_next).z+b(i,j,k_prev).z-b(i,j_prev,k_prev).z) )  -   dtdz*( by*(b(i,j,k).y-b(i,j,k_prev).y) + ay*(b(i,j_next,k).y-b(i,j_next,k_prev).y+b(i,j_prev,k).y-b(i,j_prev,k_prev).y) );
    e(i,j,k).y = e(i,j,k).y   -   dtdx*( bz*(b(i,j,k).z-b(i_prev,j,k).z) + az*(b(i,j,k_next).z-b(i_prev,j,k_next).z+b(i,j,k_prev).z-b(i_prev,j,k_prev).z) )  +   dtdz*( bx*(b(i,j,k).x-b(i,j,k_prev).x) + ax*(b(i_next,j,k).x-b(i_next,j,k_prev).x+b(i_prev,j,k).x-b(i_prev,j,k_prev).x) );
    e(i,j,k).z = e(i,j,k).z   +   dtdx*( by*(b(i,j,k).y-b(i_prev,j,k).y) + ay*(b(i,j_next,k).y-b(i_prev,j_next,k).y+b(i,j_prev,k).y-b(i_prev,j_prev,k).y) )  -   dtdy*( bx*(b(i,j,k).x-b(i,j_prev,k).x) + ax*(b(i_next,j,k).x-b(i_next,j_prev,k).x+b(i_prev,j,k).x-b(i_prev,j_prev,k).x) );
}

template <typename T>
void advance_e_periodic(T & e, const T & b) {
    const int n = e.get_n();

    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
            for(int k=0;k<n;k++) {
                int i_prev = (i != 0 ? i-1 : n-1);
                int j_prev = (j != 0 ? j-1 : n-1);
                int k_prev = (k != 0 ? k-1 : n-1);
                int i_next = (i != n-1 ? i+1 : 0);
                int j_next = (j != n-1 ? j+1 : 0);
                int k_next = (k != n-1 ? k+1 : 0);
                advance_e_point(e, b, i, j, k, i_prev, j_prev, k_prev, i_next, j_next, k_next);
            }
        }
    }
}

void advance_e_z_periodic(morton_array<vector3d> & e, const morton_array<vector3d> & b) {
    const unsigned int n = e.get_size();

    for (unsigned int curr = 0; curr < n; curr++) {
        auto x_p = b.get_x_next(curr);
        auto y_p = b.get_y_next(curr);
        auto z_p = b.get_z_next(curr);

        auto x_m = b.get_x_prev(curr);
        auto y_m = b.get_y_prev(curr);
        auto z_m = b.get_z_prev(curr);

        auto xy_mm = b.get_y_prev(x_m);
        auto xy_pm = b.get_y_prev(x_p);
        auto xy_mp = b.get_y_next(x_m);

        auto xz_mm = b.get_z_prev(x_m);
        auto xz_pm = b.get_z_prev(x_p);
        auto xz_mp = b.get_z_next(x_m);

        auto yz_mm = b.get_z_prev(y_m);
        auto yz_pm = b.get_z_prev(y_p);
        auto yz_mp = b.get_z_next(y_m);

        e[curr].x = e[curr].x   +   dtdy*( bz*(b[curr].z-b[y_m].z) + az*(b[z_p].z-b[yz_mp].z+b[z_m].z-b[yz_mm].z) )  -   dtdz*( by*(b[curr].y-b[z_m].y) + ay*(b[y_p].y-b[yz_pm].y+b[y_m].y-b[yz_mm].y) );
        e[curr].y = e[curr].y   -   dtdx*( bz*(b[curr].z-b[x_m].z) + az*(b[z_p].z-b[xz_mp].z+b[z_m].z-b[xz_mm].z) )  +   dtdz*( bx*(b[curr].x-b[z_m].x) + ax*(b[x_p].x-b[xz_pm].x+b[x_m].x-b[xz_mm].x) );
        e[curr].z = e[curr].z   +   dtdx*( by*(b[curr].y-b[x_m].y) + ay*(b[y_p].y-b[xy_mp].y+b[y_m].y-b[xy_mm].y) )  -   dtdy*( bx*(b[curr].x-b[y_m].x) + ax*(b[x_p].x-b[xy_pm].x+b[x_m].x-b[xy_mm].x) );
    }
}

template <typename T>
void advance(T & e, T & b) {
    advance_b(e, b);
    advance_e(e, b);
    advance_b(e, b);
}

void advance_z(morton_array<vector3d> & e, morton_array<vector3d> & b) {
    advance_b_z(e, b);
    advance_e_z(e, b);
    advance_b_z(e, b);
}

template <typename T>
void advance_periodic(T & e, T & b) {
    advance_b_periodic(e, b);
    advance_e_periodic(e, b);
    advance_b_periodic(e, b);
}

void advance_z_periodic(morton_array<vector3d> & e, morton_array<vector3d> & b) {
    advance_b_z_periodic(e, b);
    advance_e_z_periodic(e, b);
    advance_b_z_periodic(e, b);
}

template <typename T>
void run_test(std::function<void(T&, T&)> func, const std::string & testname, int size, int iterations=10) {
    T a(size);
    T b(size);
    vector<double> times(iterations);
    for (int i = 0; i < iterations; i++) {
        randomize(a);
        randomize(b);
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
        cout << boost::format("%15s ") % testname << boost::format("%9.3f ms") % avg
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
                if (std::fabs(b(i, j, k).x - a(i, j, k).x) > 1e-6) return false;
                if (std::fabs(b(i, j, k).y - a(i, j, k).y) > 1e-6) return false;
                if (std::fabs(b(i, j, k).z - a(i, j, k).z) > 1e-6) return false;
            }
        }
    }
    return true;
}

template <typename T>
void print_array(const T& a) {
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

void run_checks(int size) {
    cout << "Checking size " << size << "..." << endl;

    // vector3d checks

    morton_array<vector3d> e1(size);
    morton_array<vector3d> e2(size);
    morton_array<vector3d> b1(size);
    morton_array<vector3d> b2(size);

    randomize(e1);
    e2 = e1;
    randomize(b1);
    b2 = b1;

    advance_b(e1, b1);
    advance_b_z(e2, b2);

    if (!is_equal(b1, b2)) cout << "advance_b_z is different" << endl;

    randomize(e1);
    e2 = e1;
    randomize(b1);
    b2 = b1;

    advance_e(e1, b1);
    advance_e_z(e2, b2);

    if (!is_equal(e1, e2)) cout << "advance_e_z is different" << endl;

    randomize(e1);
    e2 = e1;
    randomize(b1);
    b2 = b1;

    advance(e1, b1);
    advance_z(e2, b2);

    if (!is_equal(e1, e2) || !is_equal(b1, b2)) cout << "advance_z is different" << endl;

    randomize(e1);
    e2 = e1;
    randomize(b1);
    b2 = b1;

    advance_b_periodic(e1, b1);
    advance_b_z_periodic(e2, b2);

    if (!is_equal(b1, b2)) cout << "advance_b_z_periodic is different" << endl;

    randomize(e1);
    e2 = e1;
    randomize(b1);
    b2 = b1;

    advance_e_periodic(e1, b1);
    advance_e_z_periodic(e2, b2);

    if (!is_equal(e1, e2)) cout << "advance_e_z_periodic is different" << endl;

    randomize(e1);
    e2 = e1;
    randomize(b1);
    b2 = b1;

    advance_periodic(e1, b1);
    advance_z_periodic(e2, b2);

    if (!is_equal(e1, e2) || !is_equal(b1, b2)) cout << "advance_z_periodic is different" << endl;
}

int main(int argc, char **argv) {
    dt = 0.01;
    double dx = 0.05;
    double dy = 0.05;
    double dz = 0.05;

    dtdx = dt / dx;
    dtdy = dt / dy;
    dtdz = dt / dz;

    by = 1 - kappa * dx * dx / (dy * dy);
    bz = 1 - kappa * dx * dx / (dz * dz);
    ay = 0.5 * (1 - by);
    az = 0.5 * (1 - bz);

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


        run_test<simple_array<vector3d>>(advance_b<simple_array<vector3d>>, "simple(b)", size, iterations);
        run_test<cached_array<vector3d>>(advance_b<cached_array<vector3d>>, "cached(b)", size, iterations);
        run_test<morton_array<vector3d>>(advance_b<morton_array<vector3d>>, "morton_bad(b)", size, iterations);
        run_test<morton_array<vector3d>>(advance_b_z, "morton(b)", size, iterations);

        run_test<simple_array<vector3d>>(advance_e<simple_array<vector3d>>, "simple(e)", size, iterations);
        run_test<cached_array<vector3d>>(advance_e<cached_array<vector3d>>, "cached(e)", size, iterations);
        run_test<morton_array<vector3d>>(advance_e<morton_array<vector3d>>, "morton_bad(e)", size, iterations);
        run_test<morton_array<vector3d>>(advance_e_z, "morton(e)", size, iterations);

        run_test<simple_array<vector3d>>(advance<simple_array<vector3d>>, "simple", size, iterations);
        run_test<cached_array<vector3d>>(advance<cached_array<vector3d>>, "cached", size, iterations);
        run_test<morton_array<vector3d>>(advance<morton_array<vector3d>>, "morton_bad", size, iterations);
        run_test<morton_array<vector3d>>(advance_z, "morton", size, iterations);

        run_test<simple_array<vector3d>>(advance_periodic<simple_array<vector3d>>, "simple_p", size, iterations);
        run_test<cached_array<vector3d>>(advance_periodic<cached_array<vector3d>>, "cached_p", size, iterations);
        run_test<morton_array<vector3d>>(advance_periodic<morton_array<vector3d>>, "morton_bad_p", size, iterations);
        run_test<morton_array<vector3d>>(advance_z_periodic, "morton_p", size, iterations);

    }

    MPI_Finalize();
    return 0;
}

