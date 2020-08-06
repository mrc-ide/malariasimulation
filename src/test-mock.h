/*
 * test-emergence.h
 *
 *  Created on: 5 Aug 2020
 *      Author: gc1610
 */

#ifndef SRC_TEST_MOCK_H_
#define SRC_TEST_MOCK_H_

#include <individual.h>
#include <testthat/trompeloeil.hpp>
#include "mosquito_ode.h"

class MockAPI : public ProcessAPI {
public:
    MockAPI() : ProcessAPI(
        Rcpp::XPtr<State>(static_cast<State *>(nullptr), false),
        Rcpp::XPtr<Scheduler<ProcessAPI>>(static_cast<Scheduler<ProcessAPI>*>(nullptr), false),
        Rcpp::List(),
        Rcpp::Environment()
    ) {};
    MAKE_CONST_MOCK2(get_state, const individual_index_t&(const std::string& individual, const std::string& state), override);
    MAKE_MOCK3(queue_state_update, void(const std::string& individual, const std::string& state, const std::vector<size_t>& index), override);
    MAKE_MOCK3(queue_state_update, void(const std::string& individual, const std::string& state, const individual_index_t& index), override);
    MAKE_MOCK4(queue_variable_update, void(const std::string& individual, const std::string& variable, const std::vector<size_t>& index, const variable_vector_t& value), override);
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
        0
    ) {};
    MAKE_MOCK0(get_state, state_t(), override);
};

#endif /* SRC_TEST_MOCK_H_ */
