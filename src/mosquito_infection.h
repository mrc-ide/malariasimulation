/*
 * mosquito_infection.h
 *
 *  Created on: 6 Aug 2020
 *      Author: gc1610
 */

#ifndef SRC_MOSQUITO_INFECTION_H_
#define SRC_MOSQUITO_INFECTION_H_

#include <individual.h>
#include "Random.h"

Rcpp::XPtr<process_t> create_mosquito_infection_process(
    const std::string,
    const std::string,
    const std::vector<std::string>&,
    const std::vector<std::string>&,
    const std::string&,
    RandomInterface*
);

#endif /* SRC_MOSQUITO_INFECTION_H_ */
