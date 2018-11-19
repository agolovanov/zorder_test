/*
 * array.h
 *
 *  Created on: Sep 16, 2018
 *      Author: elrond16
 */

#ifndef ARRAY_H_
#define ARRAY_H_

#include <memory>
#include <cstring>
#include "morton.h"
#include <cassert>
#include <stdexcept>

template <typename T>
class simple_array {
private:
    const int n;
    std::unique_ptr<T[]> data;

public:
    simple_array(int n) : n(n) {
        data = std::unique_ptr<T[]>(new T[n * n]);
    }

    inline T & operator()(int i, int j) const {
        return data[n * i + j];
    }

    int get_n() const { return n; }

    simple_array<T> & operator=(simple_array<T> & other) {
        memcpy(data.get(), other.data.get(), n * n * sizeof(T));

        return *this;
    }
};

template <typename T>
class cached_array {
private:
    const int n;
    std::unique_ptr<T*[]> ptrs;
    std::unique_ptr<T[]> data;

public:
    cached_array(int n) : n(n) {
        data = std::unique_ptr<T[]>(new T[n * n]);
        ptrs = std::unique_ptr<T*[]>(new T*[n]);
        for (int i = 0; i < n; i++) {
            ptrs[i] = &(data[n * i]);
        }
    }

    inline T & operator() (int i, int j) const {
        return ptrs[i][j];
    }

    int get_n() const { return n; }

    cached_array<T> & operator=(cached_array<T> & other) {
        memcpy(data.get(), other.data.get(), n * n * sizeof(T));

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

    int encode_calculate(int x) {
        int res = 0;
        for (unsigned int i = 0; i < pow; i++) {
            res += (x & (1 << i)) << i;
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
        data = std::unique_ptr<T[]>(new T[n * n]);

        cache = std::unique_ptr<unsigned int[]>(new unsigned int[n]);

        for (unsigned int i = 0; i < n; i++) {
            cache[i] = encode_calculate(i);
        }
        x_mask = cache[n-1];
        y_mask = cache[n-1] << 1;
    }

    inline T & operator() (unsigned int i, unsigned int j) const {
        return data[(cache[i] << 1) | cache[j]];
    }

    inline T & operator[] (size_t i) const {
        return data[i];
    }

    inline int get_cache(int i) const {
        return cache[i];
    }

    unsigned int get_n() const { return n; }

    unsigned int get_size() const { return n * n; }

    morton_array<T> & operator=(morton_array<T> & other) {
        memcpy(data.get(), other.data.get(), n * n * sizeof(T));

        return *this;
    }

    inline bool is_xmin(unsigned int i) const { return ((i & x_mask) == 0); }
    inline bool is_xmax(unsigned int i) const { return ((i & x_mask) == x_mask); }
    inline bool is_ymin(unsigned int i) const { return ((i & y_mask) == 0); }
    inline bool is_ymax(unsigned int i) const { return ((i & y_mask) == y_mask); }

    inline unsigned int get_x_prev(unsigned int i) const {
        return (i & y_mask) | (((i & x_mask) - 1) & x_mask);
    }

    inline unsigned int get_x_next(unsigned int i) const {
        return (i & y_mask) | (((i | y_mask) + 1) & x_mask);
    }

    inline unsigned int get_y_prev(unsigned int i) const {
        return (i & x_mask) | (((i & y_mask) - 1) & y_mask);
    }

    inline unsigned int get_y_next(unsigned int i) const {
        return (i & x_mask) | (((i | x_mask) + 1) & y_mask);
    }
};


#endif /* ARRAY_H_ */
