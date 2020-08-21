// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// All test files should include the <testthat.h>
// header file.

#include "mosquito_ode.h"
#include <testthat.h>
#include "test-mock.h"
#include "human_mortality.h"
#include "test-helpers.h"

context("Sample mothers implementation") {
  test_that("sample_mothers correctly samples mothers from the population") {
      MockRandom random;
      std::vector<bool> sampleable{true, true, false, true, true, true, true};
      std::vector<size_t> groups{1, 1, 1, 2, 2, 3, 3};
      std::vector<size_t> died{1, 3, 4};

      std::vector<size_t> first{1};
      std::vector<size_t> second{1, 0};
      std::vector<size_t> third;
      trompeloeil::sequence seq;
      REQUIRE_CALL(random, sample(2, 1))
          .IN_SEQUENCE(seq)
          .RETURN(first);
      REQUIRE_CALL(random, sample(2, 2))
          .IN_SEQUENCE(seq)
          .RETURN(second);
      REQUIRE_CALL(random, sample(2, 0))
          .IN_SEQUENCE(seq)
          .RETURN(third);

      auto mothers = sample_mothers(sampleable, groups, 3, died, &random);

      expect_true(mothers == std::vector<size_t>({1, 4, 3}));
  }

  test_that("we do not make assumptions about monotonacity") {
      MockRandom random;
      std::vector<bool> sampleable{true, true, false, true, true, true, true};
      std::vector<size_t> groups{1, 1, 1, 2, 2, 3, 3};
      std::vector<size_t> died{3, 1, 4};

      std::vector<size_t> first{1};
      std::vector<size_t> second{1, 0};
      std::vector<size_t> third;
      trompeloeil::sequence seq;
      REQUIRE_CALL(random, sample(2, 1))
          .IN_SEQUENCE(seq)
          .RETURN(first);
      REQUIRE_CALL(random, sample(2, 2))
          .IN_SEQUENCE(seq)
          .RETURN(second);
      REQUIRE_CALL(random, sample(2, 0))
          .IN_SEQUENCE(seq)
          .RETURN(third);

      auto mothers = sample_mothers(sampleable, groups, 3, died, &random);

      expect_true(mothers == std::vector<size_t>({4, 1, 3}));
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
        params["mortality_rate"] = std::vector<double>{.95};
        params["v"] = std::vector<double>{.65};
        params["fvt"] = std::vector<double>{.5};
        params["pcm"] = std::vector<double>{.774368};
        params["pvm"] = std::vector<double>{.195768};

        REQUIRE_CALL(api, get_state("human", "D"))
            .RETURN(individual_index(4, std::vector<size_t>{0}));
        REQUIRE_CALL(api, get_state("human", "Tr"))
            .RETURN(individual_index(4, std::vector<size_t>{1}));
        REQUIRE_CALL(api, get_variable("human", "is_severe"))
            .RETURN(variable_vector_t{1, 1, 0, 0});
        REQUIRE_CALL(api, get_variable("human", "zeta_group"))
            .RETURN(variable_vector_t{1, 1, 2, 2});
        REQUIRE_CALL(api, get_variable("human", "ICM"))
            .RETURN(variable_vector_t{1, 2, 3, 4});
        REQUIRE_CALL(api, get_variable("human", "IVM"))
            .RETURN(variable_vector_t{1, 2, 3, 4});
        REQUIRE_CALL(api, get_parameters())
            .RETURN(params);
        REQUIRE_CALL(api, get_timestep())
            .RETURN(2);

        // We can't easily mock sample_mothers here as we can only
        // mock out methods on classes, or write everything in a way
        // that specifically enables mocking.
        trompeloeil::sequence seq_bernoulli;
        REQUIRE_CALL(random, bernoulli(4, .95)) // mortality_rate
            .IN_SEQUENCE(seq_bernoulli)
            .RETURN(std::vector<size_t>{3});
        REQUIRE_CALL(random, bernoulli(_, .65)) // v
            .IN_SEQUENCE(seq_bernoulli)
            .RETURN(std::vector<size_t>{});
        REQUIRE_CALL(random, bernoulli(_, .5)) // fvt
            .IN_SEQUENCE(seq_bernoulli)
            .RETURN(std::vector<size_t>{0});

        trompeloeil::sequence seq_sample;
        REQUIRE_CALL(random, sample(2, 1))
            .IN_SEQUENCE(seq_sample)
            .RETURN(std::vector<size_t>{0});
        REQUIRE_CALL(random, sample(2, 1))
            .IN_SEQUENCE(seq_sample)
            .RETURN(std::vector<size_t>{3});

        std::vector<size_t> died{1, 3};
        std::vector<double> icm{.774368, 3.097472};
        std::vector<double> ivm{0.195768, 0.783072};
        std::vector<double> negone{-1};
        std::vector<double> zero{0};
        REQUIRE_CALL(api, queue_variable_update("human", "birth",
                                                died, std::vector<double>{2}));
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

        REQUIRE_CALL(api, clear_schedule("infection", died));
        REQUIRE_CALL(api, clear_schedule("asymptomatic_infection", died));

        auto mortality_process = create_mortality_process(&random);
        (*mortality_process)(api);
    }
}
