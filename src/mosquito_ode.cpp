/*
 * mosquito_ode.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "mosquito_ode.h"

integration_function_t create_ode(MosquitoModel& model) {
    return [&model](const state_t& x , state_t& dxdt , double t) {
        auto incubation_survival = exp(-model.mu * model.tau);
        dxdt[0] = model.beta * (x[3] + x[4] + x[5]) //new eggs
            - x[0] / model.de //growth to late larval stage
            - x[0] * model.mue * (1 + (x[0] + x[1]) / model.K0); //early larval deaths
        dxdt[1] = x[0] / model.de //growth from early larval
            - x[1] / model.dl //growth to pupal
            - x[1] * model.mul * (1 + model.gamma * (x[0] + x[1]) / model.K0); //late larval deaths
        dxdt[2] = x[1] / model.dl //growth to pupae
            - x[2] / model.dp //growth to adult
            - x[2] * model.mup; // death of pupae
        dxdt[3] = .5 * x[2] / model.dp //growth to adult female
            - x[3] * model.foim //infections
            - x[3] * model.mu; //deaths
        dxdt[4] = x[3] * model.foim  //infections
            - model.lagged_incubating.front() * incubation_survival //survived incubation period
            - x[4] * model.mu; // deaths
        dxdt[5] = model.lagged_incubating.front() * incubation_survival //survived incubation period
            - x[5] * model.mu; // deaths
    };
}

MosquitoModel::MosquitoModel(
    std::vector<double> init,
    double beta,
    double de,
    double mue,
    double K0,
    double gamma,
    double dl,
    double mul,
    double dp,
    double mup,
    double foim,
    double mu,
    size_t tau
    ):
    beta(beta),
    de(de),
    mue(mue),
    K0(K0),
    gamma(gamma),
    dl(dl),
    mul(mul),
    dp(dp),
    mup(mup),
    foim(foim),
    mu(mu),
    tau(tau)
    {
    for (auto i = 0u; i < tau + 1; ++i) {
        lagged_incubating.push(init[3] * foim);
    }
    auto in = state_t();
    for (auto i = 0u; i < in.size(); ++i) {
        in[i] = init[i];
    }
    ode = create_ode(*this);
    rk = boost::numeric::odeint::make_dense_output(
        a_tolerance,
        r_tolerance,
        boost::numeric::odeint::runge_kutta_dopri5<state_t>()
    );
    rk.initialize(in, 0, 1);
}

void MosquitoModel::step(double new_foim) {
    foim = new_foim;
    rk.do_step(ode);
    lagged_incubating.pop();
    lagged_incubating.push(rk.current_state()[3] * foim);
}

//[[Rcpp::export]]
Rcpp::XPtr<MosquitoModel> create_mosquito_model(
    std::vector<double> init,
    double beta,
    double de,
    double mue,
    double K0,
    double gamma,
    double dl,
    double mul,
    double dp,
    double mup,
    double foim,
    double mu,
    size_t tau
    ) {
    auto model = new MosquitoModel(
        init,
        beta,
        de,
        mue,
        K0,
        gamma,
        dl,
        mul,
        dp,
        mup,
        foim,
        mu,
        tau
    );
    return Rcpp::XPtr<MosquitoModel>(model, true);
}

//[[Rcpp::export]]
void mosquito_model_step(Rcpp::XPtr<MosquitoModel> model, double foim) {
    model->step(foim);
}

//[[Rcpp::export]]
std::vector<double> mosquito_model_get_states(Rcpp::XPtr<MosquitoModel> model) {
    const auto& state = model->rk.current_state();
    return std::vector<double>(state.cbegin(), state.cend());
}
