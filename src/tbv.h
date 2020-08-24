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
    const std::string& human,
    const std::vector<std::string>& human_states,
    const std::string& vaccinated_handle,
    const params_t& params
);

inline double calculate_tbv_antibodies(
    double t,
    double tau,
    double rho,
    double ds,
    double dl
    ) {
    return tau * (rho * exp(-t * log(2) / ds) + (1 - rho) * exp(-t * log(2) / dl));
}

inline double calculate_TRA(double mu, double gamma1, double gamma2, double antibodies) {
    auto numerator = pow(antibodies / mu, gamma1);
    return numerator / (numerator + gamma2);
}

inline double calculate_TBA(double mx, double k, double tra) {
    auto offset = pow(k / (k + mx), k);
    auto scale = 1. / (1. - offset);
    auto tra_transformation = pow(k / (k + mx * (1 - tra)), k);
    return scale * (tra_transformation - offset);
}

#endif /* SRC_TBV_H_ */
