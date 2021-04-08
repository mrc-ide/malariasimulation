/*
 * adult_mosquito_ode.h
 *
 *  Created on: 17 Mar 2021
 *      Author: gc1610
 */

#ifndef SRC_ADULT_MOSQUITO_ODE_H_
#define SRC_ADULT_MOSQUITO_ODE_H_

#include "mosquito_ode.h"
#include "history.h"

/*
 * The adult states are:
 * 3 - S  - Susceptible
 * 4 - E  - Incubating
 * 5 - I  - Infectious
 */
enum class AdultODEState : size_t {S = 3, E = 4, I = 5};

struct AdultMosquitoModel {
    MosquitoModel growth_model;
    double mu; //death rate for adult female mosquitoes
    const double tau; //extrinsic incubation period
    History susceptible; //history for susceptible counts
    History foim; //history for force of infection towards mosquitoes
    AdultMosquitoModel(MosquitoModel, double, double, double, double);
};

integration_function_t create_ode(AdultMosquitoModel& model);

#endif /* SRC_ADULT_MOSQUITO_ODE_H_ */
