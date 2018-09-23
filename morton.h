#ifndef MORTON_H_
#define MORTON_H_

#include <vector>

class morton_encoder {
private:
    int power;
    std::vector<int> cache_x;
    std::vector<int> cache_y;
    int encode_calculate(int x) {
        int res = 0;
        for (int i = 0; i < power; i++) {
            res += (x & (1 << i)) << i;
        }
        return res;
    }
public:
    morton_encoder(int power = 0) : power(power) {
        int size = 1 << power;
        cache_x = std::vector<int>(size);
        cache_y = std::vector<int>(size);

        for (int i = 0; i < size; i++) {
            cache_y[i] = encode_calculate(i);
            cache_x[i] = cache_y[i] << 1;
        }
    }

    inline int encode(int x, int y) const {
        return cache_x[x] + cache_y[y];
    }

    morton_encoder & operator=(const morton_encoder & other) {
        this->cache_x = other.cache_x;
        this->cache_y = other.cache_y;
        this->power = other.power;
        return *this;
    }
};

#endif /* MORTON_H_ */
