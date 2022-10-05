/*
 * mosquito_ode.h
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#ifndef SRC_MOSQUITO_ODE_H_
#define SRC_MOSQUITO_ODE_H_

#include <vector>
#include <cmath>
#include <type_traits>
#include "mosquito_biology.h"
#include "solver.h"

/*
 * The states are:
 * 0 - E  - Early larval stage
 * 1 - L  - Late larval stage
 * 2 - P  - Pupal stage
 */
enum class AquaticState : size_t {E, L, P};

// Provide a convenience function for getting at the index of an enumeration
template<typename T>
constexpr auto get_idx(T value)
{
    return static_cast<std::underlying_type_t<T>>(value);
}

/*
 * AquaticMosquitoModel
 *
 * A data structure for storing aquatic mosquito equation parameters for each
 * timestep
 * 
 * total_M, biting rate and mortality are updated each timestep
 */
struct AquaticMosquitoModel {
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
    const bool model_seasonality; //whether to model seasonality
    const double g0; //fourier shape parameter
    const std::vector<double> g; //fourier shape parameters
    const std::vector<double> h; //fourier shape parameters
    const double R_bar; //average rainfall
    double mum; //adult mortality rate
    double f; //biting rate
    double rainfall_floor; //minimum rainfall

    AquaticMosquitoModel(
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
        double g0,
        std::vector<double> g,
        std::vector<double> h,
        double R_bar,
        double mum,
        double f,
        double rainfall_floor
    );
    virtual ~AquaticMosquitoModel() {};
};

integration_function_t create_eqs(AquaticMosquitoModel& model);

#endif /* SRC_MOSQUITO_ODE_H_ */
