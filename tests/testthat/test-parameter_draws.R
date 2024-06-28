test_that("Out of range error", {
  p <- get_parameters()
  expect_error(set_parameter_draw(p, 0), "draw must be an integer between 1 and 1000")
  expect_error(set_parameter_draw(p, 1001), "draw must be an integer between 1 and 1000")
})

test_that("Draw overwrite works", {
  p <- get_parameters()
  p1 <- set_parameter_draw(p, 1)
  expect_identical(p1[names(parameter_draws_pf[[1]])], parameter_draws_pf[[1]])
  p2 <- set_parameter_draw(p, 1000)
  expect_identical(p2[names(parameter_draws_pf[[1000]])], parameter_draws_pf[[1000]])
})

test_that("Out of range error (vivax)", {
  p <- get_parameters(parasite = "vivax")
  expect_error(set_parameter_draw(p, 0), "draw must be an integer between 1 and 1000")
  expect_error(set_parameter_draw(p, 1001), "draw must be an integer between 1 and 1000")
})

test_that("Draw overwrite works (vivax)", {
  p <- get_parameters(parasite = "vivax")
  p1 <- set_parameter_draw(p, 1)
  expect_identical(p1[names(parameter_draws_pv[[1]])], parameter_draws_pv[[1]])
  p2 <- set_parameter_draw(p, 1000)
  expect_identical(p2[names(parameter_draws_pv[[1000]])], parameter_draws_pv[[1000]])
})
