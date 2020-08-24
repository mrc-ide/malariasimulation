/*
 * tbv.cpp
 *
 *  Created on: 20 Aug 2020
 *      Author: gc1610
 */

#include "tbv.h"

void account_for_tbv(
    variable_vector_t &infectivity,
    ProcessAPI &api,
    const std::string& human,
    const std::vector<std::string> &human_states, //{"U", "A", "D", "Tr"}
    const std::string &vaccinated_handle,
    const params_t &params) {
    auto const& vaccinated = api.get_variable(human, vaccinated_handle);
    auto mx = std::vector<double>{
        params.at("tbv_mu")[0],
        params.at("tbv_ma")[0],
        params.at("tbv_md")[0],
        params.at("tbv_mt")[0]
    };
    auto k = params.at("tbv_k")[0];
    auto tau = params.at("tbv_tau")[0];
    auto rho = params.at("tbv_rho")[0];
    auto ds = params.at("tbv_ds")[0];
    auto dl = params.at("tbv_dl")[0];
    auto mu = params.at("tbv_mu")[0];
    auto gamma1 = params.at("tbv_gamma1")[0];
    auto gamma2 = params.at("tbv_gamma2")[0];
    auto timestep = api.get_timestep();
    for (auto state_i = 0u; state_i < human_states.size(); ++state_i) {
        auto state = api.get_state(human, human_states[state_i]);
        for (auto i : state) {
            if (vaccinated[i] > -1) {
                infectivity[i] *= (1 - calculate_TBA(
                    calculate_TRA(
                        mu,
                        gamma1,
                        gamma2,
                        calculate_tbv_antibodies(timestep - vaccinated[i], tau, rho, ds, dl)
                    ),
                    mx[state_i],
                    k
                ));
            }
        }
    }
}
