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
    double foim,
    double lsm_factor
    ) : growth_model(growth_model), mu(mu), tau(tau), foim(foim), lsm_factor(lsm_factor)
{
    for (auto i = 0u; i < tau; ++i) {
        lagged_incubating.push(incubating);
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
            .5 * model.lsm_factor * x[get_idx(AquaticState::P)] / model.growth_model.dp //growth to adult female
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
    double foim,
    double lsm_factor
    ) {
    auto model = new AdultMosquitoModel(
        *growth_model,
        mu,
        tau,
        susceptible,
        foim,
        lsm_factor
    );
    return Rcpp::XPtr<AdultMosquitoModel>(model, true);
}

//[[Rcpp::export]]
void adult_mosquito_model_update(
    Rcpp::XPtr<AdultMosquitoModel> model,
    double mu,
    double foim,
    double lsm_factor,
    double susceptible,
    double f
    ) {
    model->mu = mu;
    model->foim = foim;
    model->lsm_factor = lsm_factor;
    model->growth_model.f = f;
    model->growth_model.mum = mu;
    model->lagged_incubating.push(susceptible * foim);
    if (model->lagged_incubating.size() > 0) {
        model->lagged_incubating.pop();
    }
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
