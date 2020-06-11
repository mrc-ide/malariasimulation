test_that('Simulation runs for a few timesteps', {
  sim <- run_simulation(2)
  expect_equal(nrow(sim), 2)
})
