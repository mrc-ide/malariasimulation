// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// All test files should include the <testthat.h>
// header file.

#include "mosquito_ode.h"
#include <testthat.h>
#include "test-mock.h"
#include "human_mortality.h"

context("Mortality works") {
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
