test_that("exponential_decay_process works as expected", {
  # This rate gives a halving at every timestep
  rate <- -1 / log(0.5)

  v <- individual::DoubleVariable$new(c(0,0.5,1,2,4,10))
  p <- create_exponential_decay_process(v, rate)

  individual:::execute_any_process(p, 1)
  v$.update()

  expect_equal(v$get_values(), c(0, 0.25, 0.5, 1, 2, 5))

  individual:::execute_any_process(p, 2)
  v$.update()

  expect_equal(v$get_values(), c(0, 0.125, 0.25, 0.5, 1, 2.5))
})

test_that("exponential_decay_process fails on IntegerVariable", {
  rate <- -1 / log(0.5)
  v <- individual::IntegerVariable$new(c(0,1,2,3))
  expect_error(create_exponential_decay_process(v, rate))
})
