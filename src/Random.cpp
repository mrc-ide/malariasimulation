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

std::vector<size_t> Random::sample(size_t n, size_t size, bool replacement) {
    //                                                           probs       one_based
    Rcpp::IntegerVector res = Rcpp::sample(n, size, replacement, R_NilValue, false);
    return Rcpp::as<std::vector<size_t>>(res);
}
