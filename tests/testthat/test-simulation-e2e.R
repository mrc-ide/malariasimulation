test_that('Simulation runs for a few timesteps', {
  sim <- run_simulation(100)
  expect_equal(nrow(sim), 100)
})

test_that('run_metapop_simulation fails with incorrect mixing matrix', {
  population <- 4
  timestep <- 5
  parameters <- get_parameters(list(human_population = population))
  paramlist <- list(parameters, parameters)
  # incorrect params
  mixing <- matrix(c(1, 1), nrow = 1, ncol = 2)
  expect_error(
    run_metapop_simulation(timesteps, parameters, NULL, mixing)
  )
})

test_that('run_metapop_simulation integrates two models correctly', {
  population <- 4
  timesteps <- 5
  parameters <- get_parameters(list(human_population = population))
  parametersets <- list(parameters, parameters)
  mixing <- diag(1, nrow = 2, ncol = 2)

  outputs <- run_metapop_simulation(timesteps, parametersets, NULL, mixing)
  expect_equal(length(outputs), 2)
  expect_equal(nrow(outputs[[1]]), 5)
  expect_equal(nrow(outputs[[2]]), 5)
})
