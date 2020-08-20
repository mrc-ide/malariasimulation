/*
 * tbv.cpp
 *
 *  Created on: 20 Aug 2020
 *      Author: gc1610
 */

#include "tbv.h"

void account_for_tbv(variable_vector_t &infectivity, ProcessAPI &api,
    const std::vector<std::string> &human_states,
    const std::string &vaccinated_handle, const params_t &params) {
}

void calculate_tbv_antibodies(const std::vector<double> &t,
    const params_t &params, std::vector<double> &antibodies) {
    auto tau = params.at("tbv_tau")[0];
    auto rho = params.at("tbv_rho")[0];
    auto ds = params.at("tbv_ds")[0];
    auto dl = params.at("tbv_dl")[0];
    for (auto i = 0u; i < antibodies.size(); ++i) {
        antibodies[i] = tau * (
            rho * exp(-t[i] * log(2) / ds) + (
                1 - rho
            ) * exp(-t[i] * log(2) / dl)
        );
    }
}

void calculate_TRA(const params_t &params, std::vector<double> &antibodies) {
    auto mu = params.at("tbv_tra_mu")[0];
    auto gamma1 = params.at("tbv_gamma1")[0];
    auto gamma2 = params.at("tbv_gamma2")[0];
    for (auto i = 0u; i < antibodies.size(); ++i) {
        auto numerator = pow(antibodies[i] / mu, gamma1);
        antibodies[i] = numerator / (numerator + gamma2);
    }
}

void calculate_TBA(const std::vector<double> &mx, double k,
    std::vector<double> &tra) {
    for (auto i = 0u; i < tra.size(); ++i) {
        auto power = pow(k / (k + mx[i]), k);
        auto scale = 1. / (1. - power);
        auto offset = 1. / power;
        auto tra_transformation = pow(k / (k + mx[i] * (1 - tra[i])), k);
        tra[i] = scale * (tra_transformation - offset);
    }
}
