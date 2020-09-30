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

variable_vector_t get_age(const variable_vector_t& birth, size_t t);

std::vector<double> calculate_force_of_infection(
    const variable_vector_t& age,
    const variable_vector_t& zeta,
    const variable_vector_t& infectivity,
    const params_t& parameters
);

#endif /* SRC_MOSQUITO_INFECTION_H_ */
