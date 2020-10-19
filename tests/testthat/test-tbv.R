test_that('TBV strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_tbv(
    parameters,
    start = 10,
    end = 100,
    frequency = 365,
    ages = c(1, 2, 3, 18),
    coverage = 0.8
  )
  expect_equal(parameters$tbv, TRUE)
  expect_equal(parameters$tbv_start, 10)
  expect_equal(parameters$tbv_end, 100)
  expect_equal(parameters$tbv_ages, c(1, 2, 3, 18))
  expect_equal(parameters$tbv_coverage, .8)
})
