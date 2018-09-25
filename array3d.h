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
    const int n;
    int pow;
    std::unique_ptr<T[]> data;
    std::unique_ptr<int[]> cache;

    int encode_calculate(int x) {
        int res = 0;
        for (int i = 0; i < pow; i++) {
            res += (x & (1 << i)) << (2 * i);
        }
        return res;
    }
public:
    morton_array(int n) : n(n) {
        pow = 0;
        while ((1 << pow) < n) {
            pow++;
        }
        if ((1 << pow) != n) {
            throw std::invalid_argument("The size must be a power of two");
        }
        data = std::unique_ptr<T[]>(new T[n * n * n]);

        cache = std::unique_ptr<int[]>(new int[n]);

        for (int i = 0; i < n; i++) {
            cache[i] = encode_calculate(i);
        }
    }

    inline T & operator() (int i, int j, int k) const {
        return data[(cache[i] << 2) ^ (cache[j] << 1) ^ cache[k]];
    }

    int get_n() const { return n; }

    morton_array<T> & operator=(morton_array<T> & other) {
        memcpy(data.get(), other.data.get(), n * n * n * sizeof(T));

        return *this;
    }
};



#endif /* ARRAY3D_H_ */
