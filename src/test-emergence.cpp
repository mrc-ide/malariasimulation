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
      REQUIRE_CALL(ode, get_state)
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
      std::vector<size_t>target(100);
      variable_vector_t value(100);
      for (auto i = 0ul; i < 100; ++i) {
          target[i] = i + 1000;
          value[i] = 1;
      }
      (*process)(api);
      REQUIRE_CALL(api, get_state("mosquito", "unborn"))
          .RETURN(unborn);
      REQUIRE_CALL(ode, get_state())
          .RETURN(state_t{1000, 500, 100})
          .TIMES(AT_LEAST(1));
      REQUIRE_CALL(api, queue_state_update("mosquito", "susceptible", target));
      REQUIRE_CALL(api, queue_variable_update("mosquito", "mosquito_variety", target, value));
  }

}
