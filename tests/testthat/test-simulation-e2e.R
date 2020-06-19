test_that('Simulation runs for a few timesteps', {
  sim <- run_simulation(2)
  expect_equal(nrow(sim), 2)
})

test_that('Simulation with ode runs for a few timesteps', {
  sim <- run_simulation(2, list(vector_ode = TRUE))
  expect_equal(nrow(sim), 2)
})
