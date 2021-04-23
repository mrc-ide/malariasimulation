// All test files should include the <testthat.h>
// header file.

#include "solver.h"
#include <testthat.h>
#include "test-mock.h"
#include "mosquito_emergence.h"

integration_function_t mock_integration = [](
    const state_t&,
    state_t&,
    double t
) {};

context("Emergence works") {

    test_that("emergence process fails when there are not enough individuals") {
        //mock a solver
        Solver solver(state_t{1000, 500, 100}, mock_integration);

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
            Rcpp::List::create(Rcpp::XPtr<Solver>(&solver, false)),
            Rcpp::XPtr<CategoricalVariable>(&state, false),
            Rcpp::XPtr<CategoricalVariable>(&species, false),
            {"gamb"},
            2
        );
        CATCH_REQUIRE_THROWS((*process)(1));
    }

    test_that("emergence process fails at the boundary") {
        //mock a solver
        Solver solver(state_t{1000, 500, 100}, mock_integration);

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
            Rcpp::List::create(Rcpp::XPtr<Solver>(&solver, false)),
            Rcpp::XPtr<CategoricalVariable>(&state, false),
            Rcpp::XPtr<CategoricalVariable>(&species, false),
            {"gamb"},
            2
        );
        CATCH_REQUIRE_THROWS((*process)(1));
    }

    test_that("emergence process creates the correct number of susceptibles") {
        //mock a solver
        Solver solver(state_t{1000, 500, 100}, mock_integration);

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
            Rcpp::List::create(Rcpp::XPtr<Solver>(&solver, false)),
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

        REQUIRE_CALL(state, queue_update("Sm", target));
        REQUIRE_CALL(species, queue_update("gamb", target));
        (*process)(1);
    }

    test_that("emergence process works at the boundary") {
        //mock a solver
        Solver solver(state_t{1000, 500, 100}, mock_integration);

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
            Rcpp::List::create(Rcpp::XPtr<Solver>(&solver, false)),
            Rcpp::XPtr<CategoricalVariable>(&state, false),
            Rcpp::XPtr<CategoricalVariable>(&species, false),
            {"gamb"},
            2
        );
        individual_index_t target(2000);
        for (auto i = 0ul; i < expected_to_emerge; ++i) {
            target.insert(i + 1000);
        }

        REQUIRE_CALL(state, queue_update("Sm", target));
        REQUIRE_CALL(species, queue_update("gamb", target));
        (*process)(1);
    }
}
