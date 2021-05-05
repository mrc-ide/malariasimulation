/*
 * test-history.cpp
 *
 *  Created on: 08 Apr 2021
 *      Author: gc1610
 */

#include <testthat.h>
#include "history.h"

context("History") {
    test_that("Can get history at exact time point") {
        History h;
        h.push(1., 1.);
        h.push(1.5, 4.);
        expect_true(h.at(1.) == 1.);
    }

    test_that("Can get history at interpolated time point") {
        History h;
        h.push(1., 0.);
        h.push(2., 4.);
        expect_true(h.at(1.5) == 2.);
    }

    test_that("Errors gracefully when I search before the history") {
        History h;
        h.push(1., 0.);
        h.push(2., 4.);
        expect_error(h.at(0.));
    }

    test_that("Errors gracefully when I search after the history") {
        History h;
        h.push(1., 0.);
        h.push(2., 4.);
        expect_error(h.at(4.));
    }

    test_that("Can restrict size of history") {
        History h(2);
        h.push(2.5, 3.);
        h.push(1., 0.);
        h.push(2., 4.);
        expect_error(h.at(1));
    }
}
