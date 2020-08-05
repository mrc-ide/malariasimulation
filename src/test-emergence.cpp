// All test files should include the <testthat.h>
// header file.
#include "test-mock.h"
#include "mosquito_emergence.h"

context("Emergence works") {

  test_that("emergence process fails when there are not enough individuals") {
      auto process = create_mosquito_emergence_process_cpp("human", "unborn", "susceptible");
      MockAPI api;
      (*process)(api);
  }

}
