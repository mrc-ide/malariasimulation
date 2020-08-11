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

void create_infectivity_target_vector(
    const individual_index_t&,
    const variable_vector_t&,
    const std::vector<std::vector<size_t>>&,
    std::vector<size_t>&
);

Rcpp::XPtr<process_t> create_mosquito_infection_process(
    const std::string,
    const std::string,
    const std::vector<std::string>&,
    const std::vector<std::string>&,
    RandomInterface*
);

#endif /* SRC_MOSQUITO_INFECTION_H_ */
