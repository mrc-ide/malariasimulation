test_that('Simulation runs for 10 timesteps', {
  sim <- run_simulation(2)
  expect_equal(dim(sim$states), c(100000, 2))
  expect_equal(dim(sim$variables), c(100000, 1, 2))
})
