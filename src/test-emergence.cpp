// All test files should include the <testthat.h>
// header file.

#include "mosquito_ode.h"
#include <testthat.h>
#include "test-mock.h"
#include "mosquito_emergence.h"

context("Emergence works") {

    test_that("emergence process fails when there are not enough individuals") {
        MockODE ode;
        //mock my variables
        auto state_values = std::vector<std::string>(2000);
        for (auto n = 0; n < 2000; ++n) {
            state_values[n] = "Sm";
        }
        MockCategory state(
            {"Sm", "NonExistent"},
            state_values
        );
        MockCategory species(
            {"gamb"},
            std::vector<std::string>(2000, "gamb")
        );

        auto process = create_mosquito_emergence_process_cpp(
            Rcpp::List::create(Rcpp::XPtr<MosquitoModel>(&ode, false)),
            Rcpp::XPtr<CategoricalVariable>(&state, false),
            Rcpp::XPtr<CategoricalVariable>(&species, false),
            {"gamb"},
            2
        );
        REQUIRE_CALL(ode, get_state())
        .RETURN(state_t{1000, 500, 100})
        .TIMES(AT_LEAST(1));
        CATCH_REQUIRE_THROWS((*process)(1));
    }

    test_that("emergence process fails at the boundary") {
        MockODE ode;
        //mock my variables
        auto state_values = std::vector<std::string>(2000);
        for (auto n = 0; n < 2000; ++n) {
            state_values[n] = "Sm";
        }
        for (auto n = 1000; n < 1024; ++n) {
            state_values[n] = "NonExistent";
        }
        MockCategory state(
            {"Sm", "NonExistent"},
            state_values
        );
        MockCategory species(
            {"gamb"},
            std::vector<std::string>(2000, "gamb")
        );

        auto process = create_mosquito_emergence_process_cpp(
            Rcpp::List::create(Rcpp::XPtr<MosquitoModel>(&ode, false)),
            Rcpp::XPtr<CategoricalVariable>(&state, false),
            Rcpp::XPtr<CategoricalVariable>(&species, false),
            {"gamb"},
            2
        );
        REQUIRE_CALL(ode, get_state())
        .RETURN(state_t{1000, 500, 100})
        .TIMES(AT_LEAST(1));
        CATCH_REQUIRE_THROWS((*process)(1));
    }

    test_that("emergence process creates the correct number of susceptibles") {
        MockODE ode;

        //mock my variables
        auto state_values = std::vector<std::string>(2000);
        for (auto n = 0; n < 1000; ++n) {
            state_values[n] = "Sm";
        }
        for (auto n = 1000; n < 2000; ++n) {
            state_values[n] = "NonExistent";
        }
        MockCategory state(
            {"Sm", "NonExistent"},
            state_values
        );
        MockCategory species(
            {"gamb"},
            std::vector<std::string>(2000, "gamb")
        );

        auto process = create_mosquito_emergence_process_cpp(
            Rcpp::List::create(Rcpp::XPtr<MosquitoModel>(&ode, false)),
            Rcpp::XPtr<CategoricalVariable>(&state, false),
            Rcpp::XPtr<CategoricalVariable>(&species, false),
            {"gamb"},
            2
        );
        
        auto expected_to_emerge = 25u;
        individual_index_t target(2000);
        for (auto i = 0ul; i < expected_to_emerge; ++i) {
            target.insert(i + 1000);
        }

        REQUIRE_CALL(ode, get_state())
        .RETURN(state_t{1000, 500, 100})
        .TIMES(AT_LEAST(1));
        REQUIRE_CALL(state, queue_update("Sm", target));
        REQUIRE_CALL(species, queue_update("gamb", target));
        (*process)(1);
    }

    test_that("emergence process works at the boundary") {
        MockODE ode;
        //mock my variables
        auto state_values = std::vector<std::string>(2000);
        for (auto n = 0; n < 2000; ++n) {
            state_values[n] = "Sm";
        }
        auto expected_to_emerge = 25u;
        for (auto n = 1000; n < 1000 + expected_to_emerge; ++n) {
            state_values[n] = "NonExistent";
        }
        MockCategory state(
            {"Sm", "NonExistent"},
            state_values
        );
        MockCategory species(
            {"gamb"},
            std::vector<std::string>(2000, "gamb")
        );

        auto process = create_mosquito_emergence_process_cpp(
            Rcpp::List::create(Rcpp::XPtr<MosquitoModel>(&ode, false)),
            Rcpp::XPtr<CategoricalVariable>(&state, false),
            Rcpp::XPtr<CategoricalVariable>(&species, false),
            {"gamb"},
            2
        );
        individual_index_t target(2000);
        for (auto i = 0ul; i < expected_to_emerge; ++i) {
            target.insert(i + 1000);
        }

        REQUIRE_CALL(ode, get_state())
        .RETURN(state_t{1000, 500, 100})
        .TIMES(AT_LEAST(1));
        REQUIRE_CALL(state, queue_update("Sm", target));
        REQUIRE_CALL(species, queue_update("gamb", target));
        (*process)(1);
    }
}
