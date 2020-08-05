/*
 * mosquito_emergence.cpp
 *
 *  Created on: 4 Aug 2020
 *      Author: gc1610
 */

#include "mosquito_emergence.h"

Rcpp::XPtr<process_t> create_mosquito_emergence_process_cpp(
    std::string human,
    std::string unborn,
    std::string susceptible
    ) {
    return Rcpp::XPtr<process_t>(
        new process_t([=] (ProcessAPI& api) {
            auto target_individuals = api.get_state(human, unborn);
            api.queue_state_update(human, susceptible, target_individuals);
        }),
        true
    );
}

