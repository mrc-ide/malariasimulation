/*
 * test-mosquito-infection.cpp
 *
 *  Created on: 6 Aug 2020
 *      Author: gc1610
 */

#include "mosquito_ode.h"
#include <testthat.h>
#include "tbv.h"
#include "test-mock.h"
#include "test-helper.h"

context("TBV calculations are correct") {

    test_that("Infectivity vector is modified correctly for vaccinated individuals") {
        MockAPI api;

        //Mock parameters
        params_t params;
        params["tbv_mt"] = std::vector<double>{35};
        params["tbv_md"] = std::vector<double>{46.7};
        params["tbv_ma"] = std::vector<double>{3.6};
        params["tbv_mu"] = std::vector<double>{.8};
        params["tbv_k"] = std::vector<double>{.9};

        //Mock state
        auto u = individual_index(5, std::vector<size_t>{1});
        auto a = individual_index(5, std::vector<size_t>{2});
        auto d = individual_index(5, std::vector<size_t>{3});
        auto t = individual_index(5, std::vector<size_t>{4});
        auto tbv_vaccinated = variable_vector_t{-1, -1, 50, 50, 50};

        REQUIRE_CALL(api, get_state("human", "U"))
            .RETURN(u);
        REQUIRE_CALL(api, get_state("human", "A"))
            .RETURN(a);
        REQUIRE_CALL(api, get_state("human", "D"))
            .RETURN(d);
        REQUIRE_CALL(api, get_state("human", "Tr"))
            .RETURN(t);
        REQUIRE_CALL(api, get_variable("human", "tbv_vaccinated"))
            .RETURN(tbv_vaccinated);
        REQUIRE_CALL(api, get_timestep())
            .RETURN(55);

        //execute tbv_modifier
        auto infectivity = variable_vector_t{0, .1, .15, .5, .3};
        auto human_states = std::vector<std::string>{"U", "A", "D", "Tr"};
        auto vaccinated_handle = "tbv_vaccinated";
        account_for_tbv(infectivity, api, human_states, vaccinated_handle, params);

        //check the result
        auto expected = variable_vector_t{0, .1, .0295, .351, .1930};
        expect_true(infectivity == expected);
    }

    test_that("TBV antibodies are calculated correctly") {
        params_t params;
        params["tbv_tau"] = std::vector<double>{22};
        params["tbv_rho"] = std::vector<double>{.7};
        params["tbv_ds"] = std::vector<double>{45};
        params["tbv_dl"] = std::vector<double>{591};

        auto t = std::vector<double>{0, 0, 10, 30};
        auto antibodies = std::vector<double>(4);
        calculate_tbv_antibodies(t, params, antibodies);
        auto expected = std::vector<double>{22, 22, 19.7, 16.1};
        for (auto i = 0u; i < expected.size(); ++i) {
            expect_true(antibodies[i] == Approx(expected[i]).epsilon(.1));
        }
    }

    test_that("TRAs are calculated correctly") {
        params_t params;
        params["tbv_tra_mu"] = std::vector<double>{12.63};
        params["tbv_gamma1"] = std::vector<double>{2.5};
        params["tbv_gamma2"] = std::vector<double>{.06};

        auto antibodies = std::vector<double>{22, 22, 19.7, 5};
        calculate_TRA(params, antibodies);
        auto expected = std::vector<double>{0.985, 0.985, 0.981, 0.622};
        for (auto i = 0u; i < expected.size(); ++i) {
            expect_true(antibodies[i] == Approx(expected[i]).epsilon(.001));
        }
    }

    test_that("TBAs are calculated correctly") {
        auto tra = std::vector<double>{0.685, 0.685, 0.466, 0.349};
        auto mx = std::vector<double>{35, 46.7, 3.6, .8};
        auto k = .9;

        calculate_TBA(mx, k, tra);
        auto expected = std::vector<double>{0.0638, 0.05, 0.1602, 0.2268};
        for (auto i = 0u; i < expected.size(); ++i) {
            expect_true(tra[i] == Approx(expected[i]));
        }
    }

}
