/*
 * mosquito_biology.h
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#ifndef SRC_MOSQUITO_BIOLOGY_H_
#define SRC_MOSQUITO_BIOLOGY_H_

#include <individual.h>

double carrying_capacity(
    const size_t,
    const bool,
    const double,
    const double,
    const std::vector<double>&,
    const std::vector<double>&,
    const double,
    const double
);

double rainfall(
    const size_t,
    const double,
    const double,
    const std::vector<double>&,
    const std::vector<double>&
);

double eggs_laid(double, double, double);

#endif /* SRC_MOSQUITO_BIOLOGY_H_ */
