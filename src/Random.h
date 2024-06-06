// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/*
 * Random.h
 *
 *  Created on: 6 Aug 2020
 *      Author: gc1610
 */

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

// [[Rcpp::depends(dqrng, sitmo, BH)]]
#include <vector>
#include <numeric>
#include <dqrng_generator.h>
#include <dqrng_distribution.h>

class RandomInterface {
public:
    virtual std::vector<size_t> bernoulli(size_t, double) = 0;
    virtual std::vector<size_t> bernoulli_multi_p(const std::vector<double>) = 0;
    virtual std::vector<size_t> sample(size_t, size_t, bool) = 0;
    virtual void prop_sample_bucket(
        size_t size,
        std::vector<double> probs,
        int* result
    ) = 0;
    virtual void seed(size_t) = 0;
    virtual ~RandomInterface() = default;
};

class Random : public RandomInterface {
    dqrng::rng64_t rng;
public:
    static Random& get_instance() {
        static Random instance;
        return instance;
    }
    virtual void seed(size_t);
    virtual std::vector<size_t> bernoulli(size_t, double);
    virtual std::vector<size_t> bernoulli_multi_p(const std::vector<double>);
    virtual std::vector<size_t> sample(size_t, size_t, bool);

    // Performs a weighted sample of `probs.size()` integers with:
    //  *  a distribution of `probs`
    //  *  replacement
    //  *  dqrng::xoroshiro128plus as the random number generator
    //
    // please see https://arxiv.org/pdf/1903.00227.pdf sweepingAliasTable for the
    // method.
    virtual void prop_sample_bucket(
        size_t,
        std::vector<double>,
        int*
    );
    virtual ~Random() = default;
    Random(const Random &other) = delete;
    Random(Random &&other) = delete;
    Random& operator=(const Random &other) = delete;
    Random& operator=(Random &&other) = delete;

    std::string save_state();
    void restore_state(std::string state);
private:
    Random() : rng(dqrng::generator<dqrng::xoshiro256plus>(42)) {};
};

#endif /* SRC_RANDOM_H_ */
