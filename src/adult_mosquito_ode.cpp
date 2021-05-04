/*
 * mosquito_ode.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "adult_mosquito_ode.h"

auto MAX_HISTORY_SIZE = 1000;

AdultMosquitoModel::AdultMosquitoModel(
    MosquitoModel growth_model,
    double mu,
    double tau,
    double init_susceptible,
    double init_foim
    ) :
    growth_model(growth_model),
    mu(mu),
    tau(tau),
    susceptible(History(MAX_HISTORY_SIZE)),
    foim(History(MAX_HISTORY_SIZE))
{
    //initialise some values in the history
    susceptible.push(-tau, init_susceptible);
    susceptible.push(0, init_susceptible);
    foim.push(-tau, init_foim);
    foim.push(0, init_foim);
}

observer_t create_adult_history_updater(AdultMosquitoModel& model) {
    return [&model](const state_t& x, double t) {
        model.susceptible.push(t, x[get_idx(AdultODEState::S)]);
    };
}

integration_function_t create_ode(AdultMosquitoModel& model) {
    //create original ode
    auto growth_ode = create_ode(model.growth_model);
    return [&model, growth_ode](const state_t& x, state_t& dxdt, double t) {
        //set submodel total_M
        model.growth_model.total_M =
            x[get_idx(AdultODEState::S)] +
            x[get_idx(AdultODEState::E)] +
            x[get_idx(AdultODEState::I)];

        //run the growth ode
        growth_ode(x, dxdt, t);

        //run the adult ode
        auto incubation_survival = exp(-model.mu * model.tau);
        auto lagged_incubating = model.susceptible.at(t - model.tau) * model.foim.at(t - model.tau);

        dxdt[get_idx(AdultODEState::S)] =
            .5 * x[get_idx(ODEState::P)] / model.growth_model.dp //growth to adult female
            - x[get_idx(AdultODEState::S)] * model.foim.at(t) //infections
            - x[get_idx(AdultODEState::S)] * model.mu; //deaths   

        dxdt[get_idx(AdultODEState::E)] =
            x[get_idx(AdultODEState::S)] * model.foim.at(t)  //infections
            - lagged_incubating * incubation_survival //survived incubation period
            - x[get_idx(AdultODEState::E)] * model.mu; // deaths

        dxdt[get_idx(AdultODEState::I)] = lagged_incubating * incubation_survival //survived incubation period
            - x[get_idx(AdultODEState::I)] * model.mu; // deaths
    };
}

//[[Rcpp::export]]
Rcpp::XPtr<AdultMosquitoModel> create_adult_mosquito_model(
    Rcpp::XPtr<MosquitoModel> growth_model,
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
    double f,
    size_t timestep
    ) {
    model->mu = mu;
    model->foim.push(timestep, foim);
    model->growth_model.f = f;
    model->growth_model.mum = mu;
}

//[[Rcpp::export]]
Rcpp::XPtr<Solver> create_adult_solver(
    Rcpp::XPtr<AdultMosquitoModel> model,
    std::vector<double> init) {
    return Rcpp::XPtr<Solver>(
        new Solver(
            init,
            create_ode(*model),
            create_adult_history_updater(*model)
        ),
        true
    );
}
