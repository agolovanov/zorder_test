/*
 * PerformanceTest.h
 *
 *  Created on: Sep 17, 2018
 *      Author: agolovanov
 */

#ifndef PERFORMANCETEST_H_
#define PERFORMANCETEST_H_

#include <random>
#include <functional>

template<typename T>
class PerformanceTest {
private:
    int size;
    int iterations;
    T a;
    T b;
    T c;

    std::uniform_real_distribution<double> generator { -10, 10 };
    std::default_random_engine re;

    inline double gen() {
        return generator(re);
    }

    void printArray(T&);

public:
    PerformanceTest(int size, int iterations = 10);
    void testCopyLoop(bool order);
    void randomize();
    void printArrays();
    void runTest(std::function<void(T&, T&, T&)>, std::string);
};

#endif /* PERFORMANCETEST_H_ */
