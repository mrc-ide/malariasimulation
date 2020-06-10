/*
 * mosquito_emergence.cpp
 *
 *  Created on: 15 May 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include <individual.h>

//' @title Mosquito births
//' @description
//' This is the process for mosquito birth, it defines how many new early stage
//' larvae are created on each timestep.
//' @param mosquito the mosquito individual
//' @param susceptable the susceptable mosquito state
//' @param infected the infected mosquito state
//' @param unborn the unborn mosquito state
//' @param early_larval_stage the early stage larval state
//' @param larval_growth_event the event to transition from early to late larval stage
//[[Rcpp::export]]
Rcpp::XPtr<process_t> create_egg_laying_process_cpp(
    std::string mosquito,
    std::string susceptable,
    std::string infected,
    std::string unborn,
    std::string early_larval_stage,
    std::string larval_growth_event
    ) {
    auto process = [=](ProcessAPI& api) {
        auto n_M = api.get_state(mosquito, susceptable).size() +
            api.get_state(mosquito, infected).size();
        const auto& u = api.get_state(mosquito, unborn);
        const auto& parameters = api.get_parameters();
        if (n_M > 0) {
            auto n_eggs = parameters.at("beta")[0] * n_M;
            if (n_eggs > u.size()) {
                Rcpp::stop("Run out of mosquitos");
            }
            if (n_eggs >= 1) {
                auto target = std::vector<size_t>(n_eggs);
                auto it = u.cbegin();
                for (auto i = 0; i < target.size(); ++i) {
                    target[i] = *it;
                    ++it;
                }
                api.schedule(larval_growth_event, target, parameters.at("del")[0]);
                api.queue_state_update(mosquito, early_larval_stage, target);
            }
        }
    };
    return Rcpp::XPtr<process_t>(new process_t(process));
}
