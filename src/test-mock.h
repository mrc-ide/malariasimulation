// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/*
 * test-emergence.h
 *
 *  Created on: 5 Aug 2020
 *      Author: gc1610
 */

#ifndef SRC_TEST_MOCK_H_
#define SRC_TEST_MOCK_H_

#include "mosquito_ode.h"
#include <individual.h>
#include <testthat/trompeloeil.hpp>
#include "Random.h"

class MockAPI : public ProcessAPI {
public:
    MockAPI() : ProcessAPI(
        Rcpp::XPtr<State>(static_cast<State *>(nullptr), false),
        Rcpp::XPtr<Scheduler<ProcessAPI>>(static_cast<Scheduler<ProcessAPI>*>(nullptr), false),
        Rcpp::List(),
        Rcpp::Environment()
    ) {};
    MAKE_CONST_MOCK2(get_state, const individual_index_t&(const std::string& individual, const std::string& state), override);
    MAKE_CONST_MOCK2(get_variable, const variable_vector_t&(const std::string& individual, const std::string& variable), override);
    MAKE_CONST_MOCK0(get_parameters, const params_t&(), override);
    MAKE_MOCK3(queue_state_update, void(const std::string& individual, const std::string& state, const std::vector<size_t>& index), override);
    MAKE_MOCK3(queue_state_update, void(const std::string& individual, const std::string& state, const individual_index_t& index), override);
    MAKE_MOCK4(queue_variable_update, void(const std::string& individual, const std::string& variable, const std::vector<size_t>& index, const variable_vector_t& value), override);
    MAKE_MOCK3(schedule, void(const std::string& event, const std::vector<size_t>& index, double delay), override);
    MAKE_CONST_MOCK0(get_timestep, size_t(), override);
    MAKE_MOCK2(render, void(const std::string& label, double value), override);
    MAKE_MOCK2(clear_schedule, void(const std::string& label, const std::vector<size_t>& index), override);
    MAKE_MOCK2(clear_schedule, void(const std::string& label, const individual_index_t& index), override);
};

class MockODE : public MosquitoModel {
public:
    MockODE() : MosquitoModel(
        std::vector<double>{0, 0, 0},
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        false,
        1,
        0,
        std::vector<double>{},
        std::vector<double>{},
        0
    ) {};
    MAKE_MOCK0(get_state, state_t(), override);
};

class MockRandom : public RandomInterface {
public:
    MAKE_MOCK2(bernoulli, std::vector<size_t>(size_t, double), override);
    MAKE_MOCK3(sample, std::vector<size_t>(size_t, size_t, bool), override);
};

#endif /* SRC_TEST_MOCK_H_ */
