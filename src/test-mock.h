// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/*
 * test-emergence.h
 *
 *  Created on: 5 Aug 2020
 *      Author: gc1610
 */

#ifndef SRC_TEST_MOCK_H_
#define SRC_TEST_MOCK_H_

#include "solver.h"
#include <testthat/trompeloeil.hpp>
#include <individual.h>
#include "Random.h"

struct MockCategory : public CategoricalVariable {
    MockCategory(
        const std::vector<std::string> categories,
        const std::vector<std::string> values)
        : CategoricalVariable(categories, values) {}
    MAKE_MOCK2(queue_update, void(const std::string, const individual_index_t&), override);
};

class MockRandom : public RandomInterface {
public:
    MAKE_MOCK2(bernoulli, std::vector<size_t>(size_t, double), override);
    MAKE_MOCK3(sample, std::vector<size_t>(size_t, size_t, bool), override);
};


integration_function_t mock_integration = [](
    const state_t&,
    state_t&,
    double t
) {};

#endif /* SRC_TEST_MOCK_H_ */
