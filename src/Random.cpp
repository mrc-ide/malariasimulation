// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/*
 * Random.cpp
 *
 *  Created on: 6 Aug 2020
 *      Author: gc1610
 */

#include "Random.h"
#include <Rcpp.h>

std::vector<size_t> Random::bernoulli(size_t size, double p) {
    auto successes = R::rbinom(size, p);
    Rcpp::IntegerVector indices = Rcpp::sample(size, successes, false, R_NilValue, false);
    return Rcpp::as<std::vector<size_t>>(indices);
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

size_t Random::sample_int(size_t max, std::vector<double> prob) {
    return Rcpp::sample(max, 1, false, Rcpp::wrap(prob), false)[0];
}

std::vector<double> Random::runif(size_t n) {
    return Rcpp::as<std::vector<double>>(Rcpp::runif(n));
}
