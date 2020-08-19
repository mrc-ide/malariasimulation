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
        dxdt[get_idx(ODEState::E)] = model.beta * (model.total_M) //new eggs
            - x[get_idx(ODEState::E)] / model.de //growth to late larval stage
            - x[get_idx(ODEState::E)] * model.mue * (1 + (x[get_idx(ODEState::E)] + x[get_idx(ODEState::L)]) / model.K0); //early larval deaths
        dxdt[1] = x[get_idx(ODEState::E)] / model.de //growth from early larval
            - x[get_idx(ODEState::L)] / model.dl //growth to pupal
            - x[get_idx(ODEState::L)] * model.mul * (1 + model.gamma * (x[get_idx(ODEState::E)] + x[get_idx(ODEState::L)]) / model.K0); //late larval deaths
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
    size_t total_M
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
    total_M(total_M)
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
    size_t total_M
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
        total_M
    );
    return Rcpp::XPtr<MosquitoModel>(model, true);
}

//[[Rcpp::export]]
std::vector<double> mosquito_model_get_states(Rcpp::XPtr<MosquitoModel> model) {
    auto state = model->get_state();
    return std::vector<double>(state.cbegin(), state.cend());
}

//' @title Step mosquito ODE
//' @description collects summarises the human state, sends it to the vector ode
//' and makes a step
//'
//' @param odes the models to step, one for each species
//' @param mosquito the mosquito individual handle
//' @param states a list of all of adult mosquito states (Sm, Pm, Im)
//' @param variety the handle for the mosquito variety variable
//[[Rcpp::export]]
Rcpp::XPtr<process_t> create_ode_stepping_process_cpp(
    Rcpp::List odes,
    const std::string mosquito,
    const std::vector<std::string> states,
    const std::string variety
    ) {
    return Rcpp::XPtr<process_t>(
        new process_t([=] (ProcessAPI& api) {
            const auto& v = api.get_variable(mosquito, variety);
            const auto& Sm = api.get_state(mosquito, states[0]);
            const auto& Pm = api.get_state(mosquito, states[1]);
            const auto& Im = api.get_state(mosquito, states[2]);
            auto total_M = std::vector<size_t>(odes.size(), 0);
            for (auto i = 0ull; i < v.size(); ++i) {
                if (Sm.find(i) != Sm.cend() || Pm.find(i) != Pm.cend() || Im.find(i) != Im.cend()) {
                    ++total_M[v[i] - 1];
                }
            }
            for (auto i = 0u; i < odes.size(); ++i) {
                std::stringstream label;
                label << "total_M_" << i + 1;
                api.render(label.str(), total_M[i]);
                Rcpp::as<Rcpp::XPtr<MosquitoModel>>(odes[i])->step(total_M[i]);
            }
        }),
        true
    );
}


//Exported for testing purposes
//[[Rcpp::export]]
void mosquito_model_step(Rcpp::XPtr<MosquitoModel> model, size_t total_M) {
    model->step(total_M);
}
