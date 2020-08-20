/*
 * tbv.h
 *
 *  Created on: 20 Aug 2020
 *      Author: gc1610
 */

#ifndef SRC_TBV_H_
#define SRC_TBV_H_

#include <individual.h>

void account_for_tbv(
    variable_vector_t& infectivity,
    ProcessAPI& api,
    const std::vector<std::string>& human_states,
    const std::string& vaccinated_handle,
    const params_t& params
);

void calculate_tbv_antibodies(
    const std::vector<double>& t,
    const params_t& params,
    std::vector<double>& antibodies
);

void calculate_TRA(
    const params_t& params,
    std::vector<double>& antibodies
);

void calculate_TBA(
    const std::vector<double>& mx,
    double k,
    std::vector<double>& tra
);

#endif /* SRC_TBV_H_ */
