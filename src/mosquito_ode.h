/*
 * mosquito_ode.h
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#ifndef SRC_MOSQUITO_ODE_H_
#define SRC_MOSQUITO_ODE_H_

// [[Rcpp::depends(BH)]]
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <array>
#include <queue>
#include <cmath>
#include <functional>

/*
 * The states are:
 * 0 - E  - Early larval stage
 * 1 - L  - Late larval stage
 * 2 - P  - Pupal stage
 * 3 - Sm - Susceptible female adult
 * 4 - Pm - Incubating female adult
 * 5 - Im - Infected female adult
 */
using state_t = std::array<double, 6>;
using integration_function_t = std::function<void (const state_t&, state_t&, double)>;

struct MosquitoModel {
    /* Parameters */
    const double beta; //egg laying rate
    const double de; //delay for early larval growth
    const double mue; //death rate for early larvae
    const double K0; //baseline carrying capacity
    const double gamma; //carrying capacity parameter for late larvae
    const double dl; //delay for late larval growth
    const double mul; //death rate for late larvae
    const double dp; //delay for for pupal growth
    const double mup; //death rate for pupae
    double foim; //force of infection on mosquitoes
    const double mu; //death rate for adult mosquitoes
    const size_t tau; //the delay for infection
    std::queue<double> lagged_incubating; //last tau values for incubating mosquitos

    boost::numeric::odeint::runge_kutta_dopri5<state_t> rk;
    integration_function_t ode;
    state_t inout;
    double t = 0;
    const double dt = 1;

    MosquitoModel(
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
    );
    void step(double);
};

integration_function_t create_ode(MosquitoModel& model);

#endif /* SRC_MOSQUITO_ODE_H_ */
