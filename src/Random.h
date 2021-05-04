// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/*
 * Random.h
 *
 *  Created on: 6 Aug 2020
 *      Author: gc1610
 */

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <vector>

class RandomInterface {
public:
    virtual std::vector<size_t> bernoulli(size_t, double) = 0;
    virtual std::vector<size_t> bernoulli_multi_p(const std::vector<double>) = 0;
    virtual std::vector<size_t> sample(size_t, size_t, bool) = 0;
    virtual ~RandomInterface() = default;
};

class Random : public RandomInterface {
public:
    static Random& get_instance() {
        static Random instance;
        return instance;
    }
    virtual std::vector<size_t> bernoulli(size_t, double);
    virtual std::vector<size_t> bernoulli_multi_p(const std::vector<double>);
    virtual std::vector<size_t> sample(size_t, size_t, bool);
    virtual ~Random() = default;
    Random(const Random &other) = delete;
    Random(Random &&other) = delete;
    Random& operator=(const Random &other) = delete;
    Random& operator=(Random &&other) = delete;
private:
    Random() {};
};

#endif /* SRC_RANDOM_H_ */
