#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <Rcpp.h>
#include "Random.h"

Rcpp::IntegerVector simulate_competing_outcomes(Rcpp::NumericMatrix);
Rcpp::IntegerVector simulate_competing_outcomes(RandomInterface*, Rcpp::NumericMatrix);

#endif
