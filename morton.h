/*
 * morton.h
 *
 *  Created on: Sep 21, 2018
 *      Author: agolovanov
 */

#ifndef MORTON_H_
#define MORTON_H_

int encode(int x, int y, int size) {
    int res = 0;
    for (int i = 0; i < size; i++) {
        res += (x & (1 << i)) << (i + 1);
        res += (y & (1 << i)) << i;
    }
    return res;
}

#endif /* MORTON_H_ */
