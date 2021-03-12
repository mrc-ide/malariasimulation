/*
 * mosquito_ode.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include <individual.h>
#include "mosquito_ode.h"
#include <sstream>

integration_function_t create_ode(MosquitoModel& model) {
    return [&model](const state_t& x , state_t& dxdt , double t) {
        auto K = carrying_capacity(
            t,
            model.model_seasonality,
            model.days_per_timestep,
            model.g0,
            model.g,
            model.h,
            model.K0,
            model.R_bar
        );
        dxdt[get_idx(ODEState::E)] = model.beta * (model.total_M) //new eggs
            - x[get_idx(ODEState::E)] / model.de //growth to late larval stage
            - x[get_idx(ODEState::E)] * model.mue * (1 + (x[get_idx(ODEState::E)] + x[get_idx(ODEState::L)]) / K); //early larval deaths
        dxdt[1] = x[get_idx(ODEState::E)] / model.de //growth from early larval
            - x[get_idx(ODEState::L)] / model.dl //growth to pupal
            - x[get_idx(ODEState::L)] * model.mul * (1 + model.gamma * (x[get_idx(ODEState::E)] + x[get_idx(ODEState::L)]) / K); //late larval deaths
        dxdt[2] = x[get_idx(ODEState::L)] / model.dl //growth to pupae
            - x[get_idx(ODEState::P)] / model.dp //growth to adult
            - x[get_idx(ODEState::P)] * model.mup; // death of pupae
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
    size_t total_M,
    bool model_seasonality,
    double days_per_timestep,
    double g0,
    std::vector<double> g,
    std::vector<double> h,
    double R_bar
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
    total_M(total_M),
    model_seasonality(model_seasonality),
    days_per_timestep(days_per_timestep),
    g0(g0),
    g(g),
    h(h),
    R_bar(R_bar)
    {
    for (auto i = 0u; i < state.size(); ++i) {
        state[i] = init[i];
    }
    ode = create_ode(*this);
    rk = boost::numeric::odeint::make_dense_output(
        a_tolerance,
        r_tolerance,
        boost::numeric::odeint::runge_kutta_dopri5<state_t>()
    );
}

void MosquitoModel::step(size_t new_total_M) {
    total_M = new_total_M;
    boost::numeric::odeint::integrate_adaptive(rk, ode, state, t, t + dt, dt);
    ++t;
}

state_t MosquitoModel::get_state() {
    return state;
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
    size_t total_M,
    bool model_seasonality,
    double days_per_timestep,
    double g0,
    std::vector<double> g,
    std::vector<double> h,
    double R_bar
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
        total_M,
        model_seasonality,
        days_per_timestep,
        g0,
        g,
        h,
        R_bar
    );
    return Rcpp::XPtr<MosquitoModel>(model, true);
}

//[[Rcpp::export]]
std::vector<double> mosquito_model_get_states(Rcpp::XPtr<MosquitoModel> model) {
    auto state = model->get_state();
    return std::vector<double>(state.cbegin(), state.cend());
}

//[[Rcpp::export]]
void mosquito_model_step(Rcpp::XPtr<MosquitoModel> model, size_t total_M) {
    model->step(total_M);
}
