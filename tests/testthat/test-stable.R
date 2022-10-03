test_that('simulation can run for the equilibrium', {
  population <- 1e2
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, 100)
  sim <- run_simulation_until_stable(parameters=parameters)
  expect_gt(nrow(sim), 365)
  expect_equal(mean(tail(sim$EIR_All, 365)), 100)
})

test_that('simulation exits early for max_t', {
  population <- 1e2
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, 100)
  expect_error(run_simulation_until_stable(parameters=parameters, max_t=365))
})
