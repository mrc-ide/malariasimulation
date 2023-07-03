/*
 * test-timeseries.cpp
 *
 *  Created on: 08 Apr 2021
 *      Author: gc1610
 */

#include <testthat.h>
#include "timeseries.h"

context("Timeseries") {
  // Behaviour if searching at exact provided timepoints
    test_that("Can get history at exact time point, linear = true") {
        Timeseries t;
        t.push(1., 1.);
        t.push(1.5, 4.);
        expect_true(t.at(1.) == 1.);
        expect_true(t.at(1., false) == 1.);
    }
  
    // Behaviour if searching betweeen provided timepoints
    test_that("Can get interpolated history at new time point (linear =  true)") {
        Timeseries t;
        t.push(0., 1.);
        t.push(4., 2.);
        expect_true(t.at(1.5) == 2.);
    }
    test_that("Can get last history at new time point (linear =  false)") {
      Timeseries t;
      t.push(0., 1.);
      t.push(4., 2.);
      expect_true(t.at(1.5, false) == 0.);
    }

    // Behaviour if searching before first provided timepoint
    test_that("Errors gracefully when I search before the history") {
        Timeseries t;
        t.push(0., 1.);
        t.push(4., 2.);
        expect_error(t.at(0.));
        expect_error(t.at(0., false));
    }
    
    // Behaviour if searching after last provided timepoint
    test_that("Errors gracefully when I search after the history (linear =  true)") {
        Timeseries t;
        t.push(0., 1.);
        t.push(4., 2.);
        expect_error(t.at(4.));
    }
    test_that("Can get last history when I search after the history (linear =  false)") {
      Timeseries t;
      t.push(0., 1.);
      t.push(4., 2.);
      expect_true(t.at(4., false) == 4.);
    }
  
    // Behaviour if searching outside of size
    test_that("Can restrict size of history") {
        Timeseries t(2);
        t.push(3, 2.5);
        t.push(0., 1.);
        t.push(4., 2.);
        expect_error(t.at(1));
        expect_error(t.at(1, false));
    }
}
