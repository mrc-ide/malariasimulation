/*
 * mosquito_biology.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include "mosquito_biology.h"
#include <algorithm>

//[[Rcpp::export]]
double carrying_capacity(
    const size_t timestep,
    const bool model_seasonality,
    const double days_per_timestep,
    const double g0,
    const std::vector<double>& g,
    const std::vector<double>& h,
    const double K0,
    const double R_bar
) {
    if (model_seasonality) {
        double r = rainfall(timestep, days_per_timestep, g0, g, h);
        return std::max(K0 * r / R_bar, .01);
    }
    return std::max(K0, .01);
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
    const double days_per_timestep,
    const double g0,
    const std::vector<double>& g,
    const std::vector<double>& h
) {
    double result = g0;
    for (auto i = 0u; i < g.size(); ++i) {
        result +=
            g[i] * cos(2 * M_PI * t * days_per_timestep / 365 * (i + 1)) +
            h[i] * sin(2 * M_PI * t * days_per_timestep / 365 * (i + 1));
    }
    return result;
}
