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

class Observer {
    size_t steps;
    size_t max_steps;
    integration_function_t ode;
public:
    Observer(size_t max_steps, const integration_function_t& ode) :
        steps(0u), max_steps(max_steps), ode(ode) { }
    void operator()( const state_t &x , double t )
    {
        if (++steps > max_steps) {
            Rcpp::Rcout << "steps: " << steps;
            Rcpp::Rcout << ", t: " << t;
            for (auto i = 0u; i < x.size(); ++i) {
                Rcpp::Rcout << ", x[" << i << "]:" << x[i];
            }
            state_t dx(x.size());
            ode(x, dx, t);
            for (auto i = 0u; i < dx.size(); ++i) {
                Rcpp::Rcout << ", dx[" << i << "]:" << dx[i];
            }
            Rcpp::Rcout << std::endl;
            Rcpp::stop("Too much work");
        }
    }
    void reset() {
        steps = 0;
    }
    size_t get_steps() {
        return steps;
    }
};

struct Solver {
    Solver(
        const std::vector<double>& init,
        const integration_function_t& ode,
        const double r_tol,
        const double a_tol,
        const size_t max_steps
    );
    //solver fields
    boost::numeric::odeint::dense_output_runge_kutta<
        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_dopri5<state_t>
        >
    >rk;
    void step();
    double t = 0.;
    const double dt = 1.;
    state_t state;
    integration_function_t ode;
    double r_tolerance;
    double a_tolerance;
    Observer observer;
};

#endif /* SRC_SOLVER_H_ */
