
#include "solver.h"
#include <testthat.h>
#include "test-mock.h"
#include "utils.h"

context("Competing outcomes works") {

    test_that("no outcomes") {
        Rcpp::NumericMatrix p(4, 3);
        p.fill(0);
        auto actual = simulate_competing_outcomes(p);
        auto expected = Rcpp::IntegerVector(4);
        expected.fill(NA_INTEGER);

        for (auto i = 0; i < actual.size(); ++i) {
            expect_true(actual[i] == expected[i]);
        }
    }

    test_that("some outcomes") {
        using trompeloeil::_;
        Rcpp::NumericVector data {
            0, .2, 0, .4,
            0, 0, .3, 0,
            0, 0, 0, .1
        };
        data.attr("dim") = Rcpp::Dimension(4, 3);
        Rcpp::NumericMatrix p = Rcpp::as<Rcpp::NumericMatrix>(data);
        MockRandom r;
        auto unif = std::vector<double>{0, .1, 1., .2};
        REQUIRE_CALL(r, runif(4))
            .RETURN(unif);

        trompeloeil::sequence seq;
        auto sample1 = std::vector<double>{.2, 0., 0.};
        auto sample2 = std::vector<double>{.4, 0., .1};

        REQUIRE_CALL(r, sample_int(3, sample1))
            .IN_SEQUENCE(seq)
            .RETURN(1);

        REQUIRE_CALL(r, sample_int(3, sample2))
            .IN_SEQUENCE(seq)
            .RETURN(3);

        auto actual = simulate_competing_outcomes(&r, p);
        auto expected = Rcpp::IntegerVector::create(
            NA_INTEGER,
            1,
            NA_INTEGER,
            3
        );

        for (auto i = 0; i < actual.size(); ++i) {
            expect_true(actual[i] == expected[i]);
        }
    }

}
