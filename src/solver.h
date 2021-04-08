/*
 * solver.h
 *
 *  Created on: 18 Mar 2020
 *      Author: gc1610
 */

#ifndef SRC_SOLVER_H_
#define SRC_SOLVER_H_

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
// [[Rcpp::depends(BH)]]
#include <boost/numeric/odeint.hpp>
#pragma GCC diagnostic pop

#include <vector>
#include <functional>

using state_t = std::vector<double>;
using integration_function_t = std::function<void (const state_t&, state_t&, double)>;
using observer_t = std::function<void (const state_t&, double)>;

struct Solver {
    Solver(
        const std::vector<double>& init,
        const integration_function_t& ode
    );
    Solver(
        const std::vector<double>& init,
        const integration_function_t& ode,
        const observer_t& observer
    );
    //solver fields
    boost::numeric::odeint::dense_output_runge_kutta<
        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_dopri5<state_t>
        >
    >rk;
    const double r_tolerance = 1.0e-6;
    const double a_tolerance = 1.0e-6;
    void step();
    double t = 0.;
    const double dt = 1.;
    state_t state;
    integration_function_t ode;
    observer_t observer;
};

#endif /* SRC_SOLVER_H_ */
