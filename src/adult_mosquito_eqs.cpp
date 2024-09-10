/*
 * adult_mosquito_eqs.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "adult_mosquito_eqs.h"

AdultMosquitoModel::AdultMosquitoModel(
    AquaticMosquitoModel growth_model,
    double mu,
    double tau,
    double incubating,
    double foim
    ) : growth_model(growth_model), mu(mu), tau(tau), foim(foim)
{
    for (auto i = 0u; i < tau; ++i) {
        lagged_incubating.push_back(incubating);
    }
}

integration_function_t create_eqs(AdultMosquitoModel& model) {
    //create original ode
    auto growth_eqs = create_eqs(model.growth_model);
    return [&model, growth_eqs](const state_t& x, state_t& dxdt, double t) {
        //set submodel total_M
        model.growth_model.total_M =
            x[get_idx(AdultState::S)] +
            x[get_idx(AdultState::E)] +
            x[get_idx(AdultState::I)];

        //run the growth ode
        growth_eqs(x, dxdt, t);

        //run the adult ode
        auto incubation_survival = exp(-model.mu * model.tau);

        dxdt[get_idx(AdultState::S)] =
            .5 * x[get_idx(AquaticState::P)] / model.growth_model.dp //growth to adult female
            - x[get_idx(AdultState::S)] * model.foim //infections
            - x[get_idx(AdultState::S)] * model.mu; //deaths   

        dxdt[get_idx(AdultState::E)] =
            x[get_idx(AdultState::S)] * model.foim  //infections
            - model.lagged_incubating.front() * incubation_survival //survived incubation period
            - x[get_idx(AdultState::E)] * model.mu; // deaths

        dxdt[get_idx(AdultState::I)] = model.lagged_incubating.front() * incubation_survival //survived incubation period
            - x[get_idx(AdultState::I)] * model.mu; // deaths
    };
}

//[[Rcpp::export]]
Rcpp::XPtr<AdultMosquitoModel> create_adult_mosquito_model(
    Rcpp::XPtr<AquaticMosquitoModel> growth_model,
    double mu,
    double tau,
    double susceptible,
    double foim
    ) {
    auto model = new AdultMosquitoModel(
        *growth_model,
        mu,
        tau,
        susceptible,
        foim
    );
    return Rcpp::XPtr<AdultMosquitoModel>(model, true);
}

//[[Rcpp::export]]
void adult_mosquito_model_update(
    Rcpp::XPtr<AdultMosquitoModel> model,
    double mu,
    double foim,
    double susceptible,
    double f
    ) {
    model->mu = mu;
    model->foim = foim;
    model->growth_model.f = f;
    model->growth_model.mum = mu;
    model->lagged_incubating.push_back(susceptible * foim);
    if (model->lagged_incubating.size() > 0) {
        model->lagged_incubating.pop_front();
    }
}

//[[Rcpp::export]]
std::vector<double> adult_mosquito_model_save_state(
    Rcpp::XPtr<AdultMosquitoModel> model
    ) {
    // Only the lagged_incubating needs to be saved. The rest of the model
    // state is reset at each time step by a call to update before it gets
    // used.
    return {model->lagged_incubating.begin(), model->lagged_incubating.end()};
}

//[[Rcpp::export]]
void adult_mosquito_model_restore_state(
    Rcpp::XPtr<AdultMosquitoModel> model,
    std::vector<double> state
    ) {
    model->lagged_incubating.assign(state.begin(), state.end());
}

//[[Rcpp::export]]
Rcpp::XPtr<Solver> create_adult_solver(
    Rcpp::XPtr<AdultMosquitoModel> model,
    std::vector<double> init,
    double r_tol,
    double a_tol,
    size_t max_steps
    ) {
    return Rcpp::XPtr<Solver>(
        new Solver(
            init,
            create_eqs(*model),
            r_tol,
            a_tol,
            max_steps
        ),
        true
    );
}
