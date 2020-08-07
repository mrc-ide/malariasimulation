/*
 * mosquito_infection.cpp
 *
 *  Created on: 6 Aug 2020
 *      Author: gc1610
 */

#include "mosquito_infection.h"

variable_vector_t get_age(const variable_vector_t& birth, size_t t) {
    variable_vector_t age(birth.size());
    for (auto i = 0u; i < age.size(); ++i) {
        age[i] = t - birth[i];
    }
    return age;
}

//' @title Mosquito infection process
//' @description
//' This is the process of infection for mosquitos. It results in a state
//' transition from Sm to Pm for infected mosquitos.
//'
//' NOTE: this process will become obsolete when the model is reformulated to
//' model individual mosquitos biting individual humans.
//' @param mosquito the mosquito individual
//' @param human the human individual
//' @param states a list of relevant model states (Sm, Pm)
//' @param variables a list of relevant model variables (birth, zeta, infectivity, mosquito_variety)
//' @param infection the mosquito infection event
//[[Rcpp::export]]
Rcpp::XPtr<process_t> create_mosquito_infection_process_cpp(
    const std::string mosquito,
    const std::string human,
    const std::vector<std::string>& states, //{"Sm", "Pm"}
    const std::vector<std::string>& variables, //{"birth", "zeta", "infectivity", "variety"}
    const std::string& infection
) {
    auto& random = Random::get_instance();
    return create_mosquito_infection_process(
        mosquito,
        human,
        states,
        variables,
        infection,
        &random
    );
}

Rcpp::XPtr<process_t> create_mosquito_infection_process(
    const std::string mosquito,
    const std::string human,
    const std::vector<std::string>& states, //{"Sm", "Pm"}
    const std::vector<std::string>& variables, //{"birth", "zeta", "infectivity", "variety"}
    const std::string& infection,
    RandomInterface *random
    ) {
    return Rcpp::XPtr<process_t>(
        new process_t([=] (ProcessAPI& api) {
            const auto& parameters = api.get_parameters();
            const auto age = get_age(api.get_variable(human, variables[0]), api.get_timestep());
            const auto& zeta = api.get_variable(human, variables[1]);
            const auto& infectivity = api.get_variable(human, variables[2]);
            const auto& alpha = parameters.at("blood_meal_rates");

            // Calculate force of infection per species
            auto pi = std::vector<double>(age.size());
            auto rho = parameters.at("rho")[0];
            auto a0 = parameters.at("a0")[0];
            double sum_pi = 0;
            for (auto i = 0u; i < pi.size(); ++i) {
                auto psi = 1 - rho * exp(-age[i] / a0);
                auto p = zeta[i] * psi;
                pi[i] = p;
                sum_pi += p;
            }
            double mean_infectivity = 0;
            for (auto i = 0u; i < pi.size(); ++i) {
                mean_infectivity += pi[i] * infectivity[i];
            }
            mean_infectivity /= sum_pi;
            auto lambda = std::vector<double>(alpha.size());
            for (auto i = 0u; i < lambda.size(); ++i) {
                lambda[i] = alpha[i] * mean_infectivity;
            }

            //sample mosquitos who get infected for each species
            const auto& susceptible = api.get_state(mosquito, states[0]);
            const auto& species = api.get_variable(mosquito, variables[3]);
            auto n_susceptible = std::vector<size_t>(alpha.size(), 0);
            for (auto i : susceptible) {
                ++n_susceptible[species[i] - 1];
            }

            //infected_i is the ith of each species that should become infected
            auto infected_i = std::vector<std::vector<size_t>>(alpha.size());
            auto n_infected = 0u;
            for (auto i = 0u; i < infected_i.size(); ++i) {
                auto infected = random->bernoulli(n_susceptible[i], lambda[i]);
                infected_i[i] = infected;
                n_infected += infected.size();
            }

            //create a target vector for those mosquitos
            auto target = std::vector<size_t>(n_infected);
            auto target_counter = 0u;
            auto species_counter = std::vector<size_t>(3u, 0u);
            for (auto i = 0u; i < susceptible.size(); ++i) {
                auto species_i = species[i] - 1;
                auto n_infected_species = species_counter[species_i];
                if (n_infected_species < infected_i[species_i].size() &&
                    n_infected_species == infected_i[species_i][n_infected_species]) {
                    target[target_counter] = i;
                    ++target_counter;
                }
                ++species_counter[species_i];
            }

            //set up updates
            api.queue_state_update(mosquito, states[1], target);
            api.schedule(infection, target, parameters.at("dem")[0]);
        }),
        true
    );
};
