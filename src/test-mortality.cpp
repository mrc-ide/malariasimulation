// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// All test files should include the <testthat.h>
// header file.

#include "mosquito_ode.h"
#include <testthat.h>
#include "test-mock.h"
#include "human_mortality.h"
#include "test-helper.h"

context("Sample mothers implementation") {
    test_that("sample_mothers correctly samples mothers from the population") {
        MockRandom random;
        std::vector<bool> sampleable{true, true, false, true, true, true, true};
        std::vector<double> groups{1, 1, 1, 2, 2, 3, 3};
        std::vector<size_t> died{1, 3, 4};

        std::vector<size_t> first{1};
        std::vector<size_t> second{1, 0};
        std::vector<size_t> third;
        trompeloeil::sequence seq;
        REQUIRE_CALL(random, sample(2, 1, true))
        .IN_SEQUENCE(seq)
        .RETURN(first);
        REQUIRE_CALL(random, sample(2, 2, true))
        .IN_SEQUENCE(seq)
        .RETURN(second);
        REQUIRE_CALL(random, sample(2, 0, true))
        .IN_SEQUENCE(seq)
        .RETURN(third);

        auto mothers = sample_mothers(sampleable, groups, 3, died, &random);

        expect_true(mothers == std::vector<size_t>({1, 4, 3}));
    }

    test_that("we do not make assumptions about monotonacity") {
        MockRandom random;
        std::vector<bool> sampleable{true, true, false, true, true, true, true};
        std::vector<double> groups{1, 1, 1, 2, 2, 3, 3};
        std::vector<size_t> died{3, 1, 4};

        std::vector<size_t> first{1};
        std::vector<size_t> second{1, 0};
        std::vector<size_t> third;
        trompeloeil::sequence seq;
        REQUIRE_CALL(random, sample(2, 1, true))
        .IN_SEQUENCE(seq)
        .RETURN(first);
        REQUIRE_CALL(random, sample(2, 2, true))
        .IN_SEQUENCE(seq)
        .RETURN(second);
        REQUIRE_CALL(random, sample(2, 0, true))
        .IN_SEQUENCE(seq)
        .RETURN(third);

        auto mothers = sample_mothers(sampleable, groups, 3, died, &random);

        expect_true(mothers == std::vector<size_t>({4, 1, 3}));
    }

    test_that("sample_mothers samples from the whole population if there is no one in the same het group") {
        MockRandom random;
        std::vector<bool> sampleable{true, true, false, true, false, false, false};
        std::vector<double> groups{1, 1, 1, 2, 2, 3, 3};
        std::vector<size_t> died{3, 5, 6};

        std::vector<size_t> first;
        std::vector<size_t> second{0};
        std::vector<size_t> third{1, 0};
        trompeloeil::sequence seq;
        REQUIRE_CALL(random, sample(2, 0, true))
        .IN_SEQUENCE(seq)
        .RETURN(first);
        REQUIRE_CALL(random, sample(1, 1, true))
        .IN_SEQUENCE(seq)
        .RETURN(second);
        REQUIRE_CALL(random, sample(3, 2, true))
        .IN_SEQUENCE(seq)
        .RETURN(third);

        auto mothers = sample_mothers(sampleable, groups, 3, died, &random);

        expect_true(mothers == std::vector<size_t>({3, 1, 0}));
    }

    test_that("sample_mothers takes anyone if nobody is sampleable") {
        MockRandom random;
        std::vector<bool> sampleable(7, false);
        std::vector<double> groups{1, 1, 1, 2, 2, 3, 3};
        std::vector<size_t> died{3, 5, 6};

        std::vector<size_t> s{1, 1, 0};
        REQUIRE_CALL(random, sample(7, 3, true))
        .RETURN(s);

        auto mothers = sample_mothers(sampleable, groups, 3, died, &random);

        expect_true(mothers == std::vector<size_t>({1, 1, 0}));
    }
}


context("Mortality process") {
    test_that("morality integration test") {
        using trompeloeil::_;

        // originally ported from the R

        // api <- mock_api(
        //   list(
        //     human = list(
        //       D = c(1),
        //       Tr = c(2),
        //       is_severe = c(1., 1., 0., 0.),
        //       zeta_group = c(1, 1, 2, 2),
        //       ICM = c(1, 2, 3, 4),
        //       IVM = c(1, 2, 3, 4)
        //     )
        //   ),
        //   parameters = parameters,
        //   timestep = 2
        // )

        MockAPI api;
        MockRandom random;

        // See parameters.R for values used here
        params_t params;
        params["severe_enabled"] = std::vector<double>{1};
        params["human_population"] = std::vector<double>{4};
        params["average_age"] = std::vector<double>{7000};
        params["v"] = std::vector<double>{.65};
        params["fvt"] = std::vector<double>{.5};
        params["pcm"] = std::vector<double>{.774368};
        params["pvm"] = std::vector<double>{.195768};
        params["n_heterogeneity_groups"] = std::vector<double>{2};

        auto diseased = individual_index(4, std::vector<size_t>{0});
        auto treated = individual_index(4, std::vector<size_t>{1});
        auto is_severe = variable_vector_t{1, 1, 0, 0};
        auto zeta = variable_vector_t{1, 1, 2, 2};
        auto state_ica = variable_vector_t{1, 2, 3, 4};
        auto state_iva = variable_vector_t{1, 2, 3, 4};
        auto birth = variable_vector_t{-8758, -9488, -9853, -11678};

        // NOTE: allocations should be done in the test scope
        // allocations in the mocking scope could lead to memory errors!
        REQUIRE_CALL(api, get_variable("human", "birth"))
            .RETURN(birth);
        REQUIRE_CALL(api, get_state("human", "D"))
            .RETURN(diseased);
        REQUIRE_CALL(api, get_state("human", "Tr"))
            .RETURN(treated);
        REQUIRE_CALL(api, get_variable("human", "is_severe"))
            .RETURN(is_severe);
        REQUIRE_CALL(api, get_variable("human", "zeta_group"))
            .RETURN(zeta);
        REQUIRE_CALL(api, get_variable("human", "ICA"))
            .RETURN(state_ica);
        REQUIRE_CALL(api, get_variable("human", "IVA"))
            .RETURN(state_iva);
        REQUIRE_CALL(api, get_parameters())
            .RETURN(params);
        REQUIRE_CALL(api, get_timestep())
            .RETURN(2);

        // We can't easily mock sample_mothers here as we can only
        // mock out methods on classes, or write everything in a way
        // that specifically enables mocking.
        auto three = std::vector<size_t>{3};
        std::vector<size_t> empty{};
        std::vector<size_t> zero_size_t{0};
        trompeloeil::sequence seq_bernoulli;
        REQUIRE_CALL(random, bernoulli(4, 1./7000.)) // mortality_rate
            .IN_SEQUENCE(seq_bernoulli)
            .RETURN(three);
        REQUIRE_CALL(random, bernoulli(_, .65)) // v
            .IN_SEQUENCE(seq_bernoulli)
            .RETURN(empty);
        REQUIRE_CALL(random, bernoulli(_, .5)) // fvt
            .IN_SEQUENCE(seq_bernoulli)
            .RETURN(zero_size_t);

        trompeloeil::sequence seq_sample;
        REQUIRE_CALL(random, sample(2, 1, true))
            .IN_SEQUENCE(seq_sample)
            .RETURN(zero_size_t);
        REQUIRE_CALL(random, sample(2, 1, true))
            .IN_SEQUENCE(seq_sample)
            .RETURN(zero_size_t);

        std::vector<size_t> died{3, 1};
        std::vector<double> icm{2.323104, .774368};
        std::vector<double> ivm{.587304, .195768};
        std::vector<double> negone{-1};
        std::vector<double> zero{0};
        std::vector<double> two{2};
        REQUIRE_CALL(api, queue_variable_update("human", "birth",
                                                died, two));
        REQUIRE_CALL(api, queue_variable_update("human", "last_boosted_ib",
                                                died, negone));
        REQUIRE_CALL(api, queue_variable_update("human", "last_boosted_ica",
                                                died, negone));
        REQUIRE_CALL(api, queue_variable_update("human", "last_boosted_iva",
                                                died, negone));
        REQUIRE_CALL(api, queue_variable_update("human", "last_boosted_id",
                                                died, negone));
        REQUIRE_CALL(api, queue_variable_update("human", "ICM", died, _))
            //.WITH(Approx(_2) == icm) // possibly would work? would be nice
            .WITH(Approx(_4[0]) == icm[0])
            .WITH(Approx(_4[1]) == icm[1]);
        REQUIRE_CALL(api, queue_variable_update("human", "IVM", died, _))
            .WITH(Approx(_4[0]) == ivm[0])
            .WITH(Approx(_4[1]) == ivm[1]);
        REQUIRE_CALL(api, queue_variable_update("human", "IB", died, zero));
        REQUIRE_CALL(api, queue_variable_update("human", "ICA", died, zero));
        REQUIRE_CALL(api, queue_variable_update("human", "IVA", died, zero));
        REQUIRE_CALL(api, queue_variable_update("human", "ID", died, zero));
        REQUIRE_CALL(api, queue_variable_update("human", "drug", died, zero));
        REQUIRE_CALL(api, queue_variable_update("human", "drug_time",
                                                died, negone));
        REQUIRE_CALL(api, queue_variable_update("human", "infectivity",
                                                died, zero));
        REQUIRE_CALL(api, queue_variable_update("human", "net_time",
                                                died, negone));
        REQUIRE_CALL(api, queue_variable_update("human", "spray_time",
                                                died, negone));

        REQUIRE_CALL(api, clear_schedule("infection", died));
        REQUIRE_CALL(api, clear_schedule("clinical_infection", died));
        REQUIRE_CALL(api, clear_schedule("asymptomatic_infection", died));
        REQUIRE_CALL(api, clear_schedule("subpatent_infection", died));
        REQUIRE_CALL(api, clear_schedule("recovery", died));
        REQUIRE_CALL(api, clear_schedule("throw_away_net", died));

        auto mortality_process = create_mortality_process(&random);
        (*mortality_process)(api);
    }
}
