test_that('Test model runs using falciparum', {
  parameters <- get_parameters()
  sim_res <- run_simulation(timesteps = 100, parameters = parameters)
  expect_equal(nrow(sim_res), 100)
})

test_that('Test vivax model runs', {
  parameters <- get_parameters()
  vivax_parameters <- set_parasite(parameters = parameters, parasite = "vivax")
  sim_res <- run_simulation(timesteps = 100, parameters = parameters)
  expect_equal(nrow(sim_res), 100)
})
