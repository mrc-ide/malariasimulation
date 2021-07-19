/*
 * solver.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "solver.h"

Solver::Solver(
    const std::vector<double>& init,
    const integration_function_t& ode,
    double r_tol,
    double a_tol
    ) : r_tolerance(r_tol), a_tolerance(a_tol), state(init), ode(ode)
{
    rk = boost::numeric::odeint::make_dense_output(
        a_tolerance,
        r_tolerance,
        boost::numeric::odeint::runge_kutta_dopri5<state_t>()
    );
}

void Solver::step() {
    boost::numeric::odeint::integrate_adaptive(
        rk,
        ode,
        state,
        t,
        t + dt,
        dt
    );
    ++t;
}

//[[Rcpp::export]]
std::vector<double> solver_get_states(Rcpp::XPtr<Solver> solver) {
    return solver->state;
}

//[[Rcpp::export]]
void solver_step(Rcpp::XPtr<Solver> solver) {
    solver->step();
}
