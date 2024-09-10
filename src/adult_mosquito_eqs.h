/*
 * adult_mosquito_eqs.h
 *
 *  Created on: 17 Mar 2021
 *      Author: gc1610
 */

#ifndef SRC_ADULT_MOSQUITO_EQS_H_
#define SRC_ADULT_MOSQUITO_EQS_H_

#include <queue>
#include "aquatic_mosquito_eqs.h"

/*
 * The adult states are:
 * 3 - S  - Susceptible
 * 4 - E  - Incubating
 * 5 - I  - Infectious
 */
enum class AdultState : size_t {S = 3, E = 4, I = 5};

/*
 * AdultMosquitoModel
 *
 * A data structure for storing equation parameters for each timestep
 * 
 * foim, mu and lagged_incubating are updated each timestep
 */
struct AdultMosquitoModel {
    AquaticMosquitoModel growth_model;
    std::deque<double> lagged_incubating; //last tau values for incubating mosquitos
    double mu; //death rate for adult female mosquitoes
    const double tau; //extrinsic incubation period
    double foim; //force of infection towards mosquitoes
    AdultMosquitoModel(AquaticMosquitoModel, double, double, double, double);
};

// create a system of equations for the solver
integration_function_t create_eqs(AdultMosquitoModel& model);

#endif /* SRC_ADULT_MOSQUITO_EQS_H_ */
