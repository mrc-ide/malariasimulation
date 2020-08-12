/*
 * Random.cpp
 *
 *  Created on: 6 Aug 2020
 *      Author: gc1610
 */

#include "Random.h"
#include <Rcpp.h>

std::vector<size_t> Random::bernoulli(size_t size, double p) {
    auto successes = Rcpp::rbinom(1, size, p)[0];
    Rcpp::IntegerVector indices = Rcpp::sample(size, successes, false, R_NilValue, false);
    return Rcpp::as<std::vector<size_t>>(indices);
}
