/*
 * mosquito_ode.h
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#ifndef SRC_MOSQUITO_ODE_H_
#define SRC_MOSQUITO_ODE_H_

// [[Rcpp::depends(BH)]]
#include <boost/numeric/odeint.hpp>
#include <array>
#include <queue>
#include <cmath>
#include <functional>
#include <type_traits>

/*
 * The states are:
 * 0 - E  - Early larval stage
 * 1 - L  - Late larval stage
 * 2 - P  - Pupal stage
 */
enum class ODEState : size_t {E, L, P};

// Provide a convenience function for getting at the index of an enumeration
template<typename T>
constexpr auto get_idx(T value)
{
    return static_cast<std::underlying_type_t<T>>(value);
}

using state_t = std::array<double, 3>;
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
    size_t total_M; //the number of adult female mosquitos in the model

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
        size_t total_M
    );
    virtual void step(size_t);
    virtual state_t get_state();
    virtual ~MosquitoModel() {};
private:
    //solver fields
    boost::numeric::odeint::dense_output_runge_kutta<
        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_dopri5<state_t>
        >
    >rk;
    const double r_tolerance = 1.0e-6;
    const double a_tolerance = 1.0e-6;
    integration_function_t ode;

    double t = 0.;
    const double dt = 1.;
    state_t state;
};

integration_function_t create_ode(MosquitoModel& model);

#endif /* SRC_MOSQUITO_ODE_H_ */
