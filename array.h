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

template <typename T>
class simple_array {
private:
    const int n1;
    const int n2;
    std::unique_ptr<T[]> data;

public:
    simple_array(int n1, int n2) : n1(n1), n2(n2) {
        data = std::unique_ptr<T[]>(new T[n1 * n2]);
    }

    inline T & operator()(int i, int j) const {
        return data[n2 * i + j];
    }

    int get_n1() const { return n1; }
    int get_n2() const { return n2; }

    simple_array<T> & operator=(simple_array<T> & other) {
        memcpy(data.get(), other.data.get(), n1 * n2 * sizeof(T));

        return *this;
    }
};

template <typename T>
class cached_array {
private:
    const int n1;
    const int n2;
    std::unique_ptr<T*[]> ptrs;
    std::unique_ptr<T[]> data;

public:
    cached_array(int n1, int n2) : n1(n1), n2(n2) {
        data = std::unique_ptr<T[]>(new T[n1 * n2]);
        ptrs = std::unique_ptr<T*[]>(new T*[n1]);
        for (int i = 0; i < n1; i++) {
            ptrs[i] = &(data[n2 * i]);
        }
    }

    inline T & operator() (int i, int j) const {
        return ptrs[i][j];
    }

    int get_n1() const { return n1; }
    int get_n2() const { return n2; }

    cached_array<T> & operator=(cached_array<T> & other) {
        memcpy(data.get(), other.data.get(), n1 * n2 * sizeof(T));

        return *this;
    }
};



#endif /* ARRAY_H_ */
