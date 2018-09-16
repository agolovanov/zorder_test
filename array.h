/*
 * array.h
 *
 *  Created on: Sep 16, 2018
 *      Author: elrond16
 */

#ifndef ARRAY_H_
#define ARRAY_H_

#include <memory>

template <typename T>
class simple_array {
private:
    int n1;
    int n2;
    std::unique_ptr<T[]> data;

public:
    simple_array(int n1, int n2) : n1(n1), n2(n2) {
        data = std::unique_ptr<T[]>(new T[n1 * n2]);
    }

    inline T & operator() (int i, int j) {
        return data[n1 * i + j];
    }
};

template <typename T>
class cached_array {
private:
    int n1;
    int n2;
    std::unique_ptr<T*[]> ptrs;
    std::unique_ptr<T[]> data;

public:
    cached_array(int n1, int n2) : n1(n1), n2(n2) {
        data = std::unique_ptr<T[]>(new T[n1 * n2]);
        ptrs = std::unique_ptr<T*[]>(new T*[n1]);
        for (int i = 0; i < n1; i++) {
            ptrs[i] = &(data[n1 * i]);
        }
    }

    inline T & operator() (int i, int j) {
        return ptrs[i][j];
    }
};



#endif /* ARRAY_H_ */
