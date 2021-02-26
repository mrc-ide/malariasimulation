test_that('Simulation runs for a few timesteps', {
  sim <- run_simulation(100)
  expect_equal(nrow(sim), 100)
})
