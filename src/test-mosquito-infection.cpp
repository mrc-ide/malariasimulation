/*
 * test-mosquito-infection.cpp
 *
 *  Created on: 6 Aug 2020
 *      Author: gc1610
 */

#include "mosquito_ode.h"
#include <testthat.h>
#include "mosquito_infection.h"
#include "test-mock.h"
#include <individual.h>

individual_index_t individual_index(size_t size, std::vector<size_t> i) {
    return individual_index_t(size, i.cbegin(), i.cend());
}

context("Mosquito infection works") {

    test_that("Mosquito infection creates the correct updates") {
        using trompeloeil::_;

        MockRandom random;
        MockAPI api;
        auto process = create_mosquito_infection_process(
            "mosquito",
            "human",
            std::vector<std::string>{"Sm", "Pm"},
            std::vector<std::string>{"birth", "zeta", "infectivity", "mosquito_variety"},
            "mosquito_infection",
            &random
        );

        //Mock state
        auto sm = individual_index(4, std::vector<size_t>{1, 2, 3, 4});
        auto birth = variable_vector_t{-7295,  -8755,  -1820, -14230};
        auto v = variable_vector_t{1, 2, 3, 3};
        auto zeta = variable_vector_t{.2, .3, .5, .9};
        auto infectivity = variable_vector_t{.2, 0, .5, .3};
        REQUIRE_CALL(api, get_variable("human", "birth"))
            .RETURN(birth);
        REQUIRE_CALL(api, get_state("mosquito", "Sm"))
            .RETURN(sm);
        REQUIRE_CALL(api, get_variable("mosquito", "mosquito_variety"))
            .RETURN(v);
        REQUIRE_CALL(api, get_variable("human", "zeta"))
            .RETURN(zeta);
        REQUIRE_CALL(api, get_variable("human", "infectivity"))
            .RETURN(infectivity);
        REQUIRE_CALL(api, get_timestep())
            .RETURN(5);

        //Mock parameters
        params_t params;
        params["blood_meal_rates"] = std::vector<double>{.92, .74, .94};
        params["dem"] = std::vector<double>{5};
        params["rho"] = std::vector<double>{.85};
        params["a0"] = std::vector<double>{8 * 365};
        REQUIRE_CALL(api, get_parameters())
            .RETURN(params);

        //Mock random
        //Check that correct lambdas are calculated for each species
        trompeloeil::sequence seq;
        auto none = std::vector<size_t>();
        auto first = std::vector<size_t>{0};
        auto both = std::vector<size_t>{0, 1};
        REQUIRE_CALL(random, bernoulli(0, _))
            .IN_SEQUENCE(seq)
            .WITH(Approx(_2) == 0.2477872)
            .RETURN(none);
        REQUIRE_CALL(random, bernoulli(1, _))
            .IN_SEQUENCE(seq)
            .WITH(Approx(_2) == 0.1993071)
            .RETURN(first);
        REQUIRE_CALL(random, bernoulli(2, _))
            .IN_SEQUENCE(seq)
            .WITH(Approx(_2) == 0.2531739)
            .RETURN(both);

        //Mock updates
        //Check that correct updates are made and scheduled
        auto state_update = std::vector<size_t>{1, 2, 3};
        REQUIRE_CALL(api, queue_state_update("mosquito", "Pm", state_update));
        REQUIRE_CALL(api, schedule("mosquito_infection", state_update, 5));

        (*process)(api);
    }

}
