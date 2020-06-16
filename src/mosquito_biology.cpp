/*
 * mosquito_biology.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include "mosquito_biology.h"

// Seasonality not yet supported
double carrying_capacity(const size_t timestep, const params_t& parameters) {
    return parameters.at("K0")[0];
}
