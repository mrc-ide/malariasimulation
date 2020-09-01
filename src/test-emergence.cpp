// All test files should include the <testthat.h>
// header file.

#include "mosquito_ode.h"
#include <testthat.h>
#include "test-mock.h"
#include "mosquito_emergence.h"

context("Emergence works") {

    test_that("emergence process fails when there are not enough individuals") {
        MockODE ode;
        auto process = create_mosquito_emergence_process_cpp(
            "mosquito",
            Rcpp::List::create(Rcpp::XPtr<MosquitoModel>(&ode, false)),
            "unborn",
            "susceptible",
            "mosquito_variety",
            2
        );
        MockAPI api;
        const individual_index_t unborn(2000);
        REQUIRE_CALL(api, get_state("mosquito", "unborn"))
        .RETURN(unborn);
        REQUIRE_CALL(ode, get_state())
        .RETURN(state_t{1000, 500, 100})
        .TIMES(AT_LEAST(1));
        CATCH_REQUIRE_THROWS((*process)(api));
    }

    test_that("emergence process fails at the boundary") {
        MockODE ode;
        auto process = create_mosquito_emergence_process_cpp(
            "mosquito",
            Rcpp::List::create(Rcpp::XPtr<MosquitoModel>(&ode, false)),
            "unborn",
            "susceptible",
            "mosquito_variety",
            2
        );
        MockAPI api;
        individual_index_t unborn(2000);
        for (auto n = 1000ul; n < 1024; ++n) {
            unborn.insert(n);
        }
        REQUIRE_CALL(api, get_state("mosquito", "unborn"))
        .RETURN(unborn);
        REQUIRE_CALL(ode, get_state())
        .RETURN(state_t{1000, 500, 100})
        .TIMES(AT_LEAST(1));
        CATCH_REQUIRE_THROWS((*process)(api));
    }

    test_that("emergence process creates the correct number of susceptibles") {
        MockODE ode;
        auto process = create_mosquito_emergence_process_cpp(
            "mosquito",
            Rcpp::List::create(Rcpp::XPtr<MosquitoModel>(&ode, false)),
            "unborn",
            "susceptible",
            "mosquito_variety",
            2
        );
        MockAPI api;
        individual_index_t unborn(2000);
        for (auto n = 1000ul; n < 2000; ++n) {
            unborn.insert(n);
        }
        auto expected_to_emerge = 25u;
        std::vector<size_t>target(expected_to_emerge);
        variable_vector_t value(expected_to_emerge);
        for (auto i = 0ul; i < expected_to_emerge; ++i) {
            target[i] = i + 1000;
            value[i] = 1;
        }
        REQUIRE_CALL(api, get_state("mosquito", "unborn"))
        .RETURN(unborn);
        REQUIRE_CALL(ode, get_state())
        .RETURN(state_t{1000, 500, 100})
        .TIMES(AT_LEAST(1));
        REQUIRE_CALL(api, queue_state_update("mosquito", "susceptible", target));
        REQUIRE_CALL(api, queue_variable_update("mosquito", "mosquito_variety", target, value));
        (*process)(api);
    }

    test_that("emergence process works at the boundary") {
        MockODE ode;
        auto process = create_mosquito_emergence_process_cpp(
            "mosquito",
            Rcpp::List::create(Rcpp::XPtr<MosquitoModel>(&ode, false)),
            "unborn",
            "susceptible",
            "mosquito_variety",
            2
        );
        MockAPI api;
        individual_index_t unborn(2000);
        auto expected_to_emerge = 25u;
        for (auto n = 1000ul; n < 1000 + expected_to_emerge; ++n) {
            unborn.insert(n);
        }
        std::vector<size_t>target(expected_to_emerge);
        variable_vector_t value(expected_to_emerge);
        for (auto i = 0ul; i < expected_to_emerge; ++i) {
            target[i] = i + 1000;
            value[i] = 1;
        }
        REQUIRE_CALL(api, get_state("mosquito", "unborn"))
        .RETURN(unborn);
        REQUIRE_CALL(ode, get_state())
        .RETURN(state_t{1000, 500, 100})
        .TIMES(AT_LEAST(1));
        REQUIRE_CALL(api, queue_state_update("mosquito", "susceptible", target));
        REQUIRE_CALL(api, queue_variable_update("mosquito", "mosquito_variety", target, value));
        (*process)(api);
    }
}
