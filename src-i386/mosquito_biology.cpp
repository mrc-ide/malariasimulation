/*
 * mosquito_biology.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <cmath>
#include "mosquito_biology.h"

//[[Rcpp::export]]
double carrying_capacity(
    const size_t timestep,
    const bool model_seasonality,
    const double g0,
    const std::vector<double>& g,
    const std::vector<double>& h,
    const double K0,
    const double R_bar,
    const double rainfall_floor
) {
    if (model_seasonality) {
        double r = rainfall(timestep, g0, g, h, rainfall_floor);
        return K0 * r / R_bar;
    }
    return K0;
}

//[[Rcpp::export]]
double eggs_laid(
    double beta,
    double mu,
    double f
) {
    auto eov = beta / mu * (exp(mu / f) - 1);
    return eov * mu * exp(-mu / f) / (1 - exp(-mu / f));
}

//[[Rcpp::export]]
double rainfall(
    const size_t t,
    const double g0,
    const std::vector<double>& g,
    const std::vector<double>& h,
    const double floor
) {
    double result = g0;
    for (auto i = 0u; i < g.size(); ++i) {
        result +=
            g[i] * cos(2 * M_PI * t * (i + 1) / 365) +
            h[i] * sin(2 * M_PI * t * (i + 1) / 365);
    }
    return std::max(result, floor);
}
