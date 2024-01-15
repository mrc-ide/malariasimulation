// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/*
 * Random.cpp
 *
 *  Created on: 6 Aug 2020
 *      Author: gc1610
 */

#include "Random.h"
#include <Rcpp.h>

void Random::seed(size_t seed) {
    rng = dqrng::generator<dqrng::xoshiro256plus>(seed);
}

std::vector<size_t> Random::bernoulli(size_t size, double p) {
    auto successes = R::rbinom(size, p);
    return sample(size, successes, false);
}

std::vector<size_t> Random::bernoulli_multi_p(const std::vector<double> p) {
    auto successes = std::vector<size_t>();
    Rcpp::NumericVector r = Rcpp::runif(p.size());
    for (auto i = 0u; i < p.size(); ++i) {
        if (r[i] < p[i]) {
            successes.push_back(i);
        }
    }
    return successes;
}

std::vector<size_t> Random::sample(size_t n, size_t size, bool replacement) {
    //                                                           probs       one_based
    Rcpp::IntegerVector res = Rcpp::sample(n, size, replacement, R_NilValue, false);
    return Rcpp::as<std::vector<size_t>>(res);
}

void Random::prop_sample_bucket(
        size_t size,
        std::vector<double> probs,
        int* result
    ) {
    auto n = probs.size();
    auto total = std::accumulate(
        probs.begin(),
        probs.end(),
        0.
    );

    // create alias table
    auto dividing_probs = probs;
    auto alternative_index = std::vector<size_t>(n);

    auto bucket_weight = total / n;

    // get first heavy
    auto heavy = 0u;
    while (heavy < n && probs[heavy] <= bucket_weight)
        ++heavy;

    // all probabilities are the same
    if (heavy == n) {
        for (auto i = 0; i < size; ++i) {
            *result = (*rng)(n);
            ++result;
        }
        return;
    }

    // get first light
    auto light = 0u;
    while (probs[light] > bucket_weight)
        ++light;

    auto residual = probs[heavy];
    size_t next_heavy;
    double packing_weight;
    while (true) {
        if (residual > bucket_weight) {
            // check if we're done
            if (light == n) {
                break;
            }

            // pack a light bucket with the residual
            alternative_index[light] = heavy;

            // update residual
            packing_weight = (light == n) ? 0 : probs[light];
            residual = residual + packing_weight - bucket_weight;

            // find the next light element
            ++light;
            while(light < n && probs[light] > bucket_weight)
                ++light;
        } else {
            // find the next heavy
            next_heavy = heavy + 1;
            while(next_heavy < n && probs[next_heavy] <= bucket_weight)
                ++next_heavy;

            dividing_probs[heavy] = residual;

            // check if we're done
            if (next_heavy == n) {
                alternative_index[heavy] = heavy;
                break;
            }

            // pack this (ex-)heavy bucket with the next heavy
            alternative_index[heavy] = next_heavy;
            // update the residual for the next heavy
            residual = residual + probs[next_heavy] - bucket_weight;
            heavy = next_heavy;
        }
    }

    // normalise the dividing_probs
    for (auto& p : dividing_probs) {
        p /= bucket_weight;
    }

    // sample
    for (auto i = 0; i < size; ++i) {
        size_t bucket = (*rng)(n);
        double acceptance = dqrng::uniform01((*rng)());
        *result = (acceptance < dividing_probs[bucket]) ? bucket :
            alternative_index[bucket];
        ++result;
    }
}

std::string Random::save_state() {
    std::ostringstream stream;
    stream << *rng;
    return stream.str();
}

void Random::restore_state(std::string state) {
    std::istringstream stream(state);
    stream >> *rng;
}
