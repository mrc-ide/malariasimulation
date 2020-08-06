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
        dxdt[0] = model.beta * (model.total_M) //new eggs
            - x[0] / model.de //growth to late larval stage
            - x[0] * model.mue * (1 + (x[0] + x[1]) / model.K0); //early larval deaths
        dxdt[1] = x[0] / model.de //growth from early larval
            - x[1] / model.dl //growth to pupal
            - x[1] * model.mul * (1 + model.gamma * (x[0] + x[1]) / model.K0); //late larval deaths
        dxdt[2] = x[1] / model.dl //growth to pupae
            - x[2] / model.dp //growth to adult
            - x[2] * model.mup; // death of pupae
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

void MosquitoModel::step(size_t new_total_M) {
    total_M = new_total_M;
    rk.do_step(ode);
}

state_t MosquitoModel::get_state() {
    return rk.current_state();
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
void mosquito_model_step(Rcpp::XPtr<MosquitoModel> model, size_t total_M) {
    model->step(total_M);
}

//[[Rcpp::export]]
std::vector<double> mosquito_model_get_states(Rcpp::XPtr<MosquitoModel> model) {
    auto state = model->get_state();
    return std::vector<double>(state.cbegin(), state.cend());
}
