/*
 * array3d.h
 *
 *  Created on: Sep 23, 2018
 *      Author: elrond16
 */

#ifndef ARRAY3D_H_
#define ARRAY3D_H_

#include <memory>
#include <cstring>
#include <cassert>
#include <stdexcept>

template <typename T>
class simple_array {
private:
    const int n;
    std::unique_ptr<T[]> data;

public:
    simple_array(int n) : n(n) {
        data = std::unique_ptr<T[]>(new T[n * n * n]);
    }

    inline T & operator()(int i, int j, int k) const {
        return data[n * n * i + n * j + k];
    }

    int get_n() const { return n; }

    simple_array<T> & operator=(simple_array<T> & other) {
        memcpy(data.get(), other.data.get(), n * n * n * sizeof(T));

        return *this;
    }
};

template <typename T>
class cached_array {
private:
    const int n;
    std::unique_ptr<T*[]> ptrs_inner;
    std::unique_ptr<T**[]> ptrs;
    std::unique_ptr<T[]> data;

public:
    cached_array(int n) : n(n) {
        data = std::unique_ptr<T[]>(new T[n * n * n]);
        ptrs_inner = std::unique_ptr<T*[]>(new T*[n * n]);
        ptrs = std::unique_ptr<T**[]>(new T**[n]);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                ptrs_inner[n * i + j] = &(data[n * n * i + n * j]);
            }
        }
        for (int i = 0; i < n; i++) {
            ptrs[i] = &(ptrs_inner[n * i]);
        }
    }

    inline T & operator() (int i, int j, int k) const {
        return ptrs[i][j][k];
    }

    int get_n() const { return n; }

    cached_array<T> & operator=(cached_array<T> & other) {
        memcpy(data.get(), other.data.get(), n * n * n * sizeof(T));

        return *this;
    }
};

template <typename T>
class morton_array {
private:
    const unsigned int n;
    unsigned int pow;
    std::unique_ptr<T[]> data;
    std::unique_ptr<unsigned int[]> cache;
    unsigned int x_mask;
    unsigned int y_mask;
    unsigned int z_mask;
    unsigned int xy_mask;
    unsigned int yz_mask;
    unsigned int xz_mask;

    int encode_calculate(unsigned int x) {
        unsigned int res = 0;
        for (unsigned int i = 0; i < pow; i++) {
            res += (x & (1 << i)) << (2 * i);
        }
        return res;
    }
public:
    morton_array(unsigned int n) : n(n) {
        pow = 0;
        while ((1u << pow) < n) {
            pow++;
        }
        if ((1u << pow) != n) {
            throw std::invalid_argument("The size must be a power of two");
        }
        data = std::unique_ptr<T[]>(new T[n * n * n]);

        cache = std::unique_ptr<unsigned int[]>(new unsigned int[n]);

        for (unsigned int i = 0; i < n; i++) {
            cache[i] = encode_calculate(i);
        }
        x_mask = cache[n-1];
        y_mask = cache[n-1] << 1;
        z_mask = cache[n-1] << 2;
        xy_mask = x_mask | y_mask;
        xz_mask = x_mask | z_mask;
        yz_mask = y_mask | z_mask;
    }

    inline T & operator() (unsigned int i, unsigned int j, unsigned int k) const {
        return data[(cache[i] << 2) | (cache[j] << 1) | cache[k]];
    }

    inline T & operator[] (size_t i) const {
        return data[i];
    }

    unsigned int get_n() const { return n; }

    unsigned int get_size() const { return n * n * n; }

    morton_array<T> & operator=(morton_array<T> & other) {
        memcpy(data.get(), other.data.get(), get_size() * sizeof(T));

        return *this;
    }

    inline bool is_xmin(unsigned int i) const { return ((i & x_mask) == 0); }
    inline bool is_xmax(unsigned int i) const { return ((i & x_mask) == x_mask); }
    inline bool is_ymin(unsigned int i) const { return ((i & y_mask) == 0); }
    inline bool is_ymax(unsigned int i) const { return ((i & y_mask) == y_mask); }
    inline bool is_zmin(unsigned int i) const { return ((i & z_mask) == 0); }
    inline bool is_zmax(unsigned int i) const { return ((i & z_mask) == z_mask); }

    inline int get_x_prev(unsigned int i) const {
        return (i & yz_mask) | (((i & x_mask) - 1) & x_mask);
    }

    inline int get_x_next(unsigned int i) const {
        return (i & yz_mask) | (((i | yz_mask) + 1) & x_mask);
    }

    inline int get_y_prev(unsigned int i) const {
        return (i & xz_mask) | (((i & y_mask) - 1) & y_mask);
    }

    inline int get_y_next(unsigned int i) const {
        return (i & xz_mask) | (((i | xz_mask) + 1) & y_mask);
    }

    inline int get_z_prev(unsigned int i) const {
        return (i & xy_mask) | (((i & z_mask) - 1) & z_mask);
    }

    inline int get_z_next(unsigned int i) const {
        return (i & xy_mask) | (((i | xy_mask) + 1) & z_mask);
    }
};



#endif /* ARRAY3D_H_ */
