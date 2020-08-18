// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/*
 * test-random.cpp
 *
 *  Created on: 12 Aug 2020
 *      Author: gc1610
 */


#include <testthat.h>
#include <Rcpp.h>
#include <unordered_set>
#include "Random.h"

void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}

context("Random is sane") {
    test_that("Bernoulli produces the right size, range of values, no duplicates") {
        auto& random = Random::get_instance();
        set_seed(42);
        auto indices = random.bernoulli(1000, .5);
        expect_true(indices.size() == 486u);
        for (auto i : indices) {
            expect_true(i < 1000u);
            expect_true(i >= 0u);
        }
        auto values = std::unordered_set<size_t>(indices.begin(), indices.end());
        expect_true(indices.size() == values.size());
    }

}

context("Sample is sane") {
    test_that("Sample produces the right size") {
        auto& random = Random::get_instance();
        set_seed(42);
        auto indices = random.sample(1000, 500);
        expect_true(indices.size() == 500u);
        for (auto i : indices) {
            expect_true(i < 1000u);
            expect_true(i >= 0u);
        }
        auto values = std::unordered_set<size_t>(indices.begin(), indices.end());
        expect_true(indices.size() == values.size());
    }
}
