/*
 * mosquito_emergence.cpp
 *
 *  Created on: 4 Aug 2020
 *      Author: gc1610
 */

#include "mosquito_emergence.h"
#include "mosquito_ode.h"

//' @title Mosquito emergence process
//' @description Move mosquitos from Unborn to Sm in line with the number of
//' pupals in the ODE models
//'
//' @param mosquito the handle for the mosquito individual
//' @param odes a list of odes for each species
//' @param unborn the handle for the unborn mosquito state
//' @param susceptible the handle for the susceptible mosquito state
//' @param variety the handle for the variable representing the mosquito species
//' @param dpl the delay for pupal growth (in timesteps)
//[[Rcpp::export]]
Rcpp::XPtr<process_t> create_mosquito_emergence_process_cpp(
    std::string mosquito,
    Rcpp::List odes,
    std::string unborn,
    std::string susceptible,
    std::string variety,
    double dpl
    ) {
    auto rate = .5 * 1./dpl;
    return Rcpp::XPtr<process_t>(
        new process_t([=] (ProcessAPI& api) {
            auto n = 0u;
            for (Rcpp::XPtr<MosquitoModel> ode : odes) {
                n += ode->get_state()[2] * rate;
            }
            auto source = api.get_state(mosquito, unborn);
            if (source.size() < n) {
                Rcpp::stop("Not enough mosquitos. Please raise mosquito_limit");
            }
            if (n > 0) {
                variable_vector_t species(n);
                auto ode_i = 1;
                auto species_i = 0u;
                for (Rcpp::XPtr<MosquitoModel> ode : odes) {
                    auto to_set = static_cast<size_t>(ode->get_state()[2] * rate);
                    auto start = species_i;
                    for (;species_i < start + to_set; ++species_i) {
                        species[species_i] = ode_i;
                    }
                    ++ode_i;
                }
                std::vector<size_t> target(n);
                auto sourceit = source.begin();
                for (auto i = 0; i < target.size(); ++i) {
                    target[i] = *sourceit;
                    ++sourceit;
                }
                api.queue_state_update(mosquito, susceptible, target);
                api.queue_variable_update(mosquito, variety, target, species);
            }
        }),
        true
    );
}
