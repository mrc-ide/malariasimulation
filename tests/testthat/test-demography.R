test_that('calculate_initial_ages defaults to an exponential distribution', {
  parameters <- get_parameters(list(human_population = 4))
  mock_exp <- mockery::mock(c(0, 1, 2, 3))
  mockery::stub(calculate_initial_ages, 'rexp', mock_exp)
  ages <- calculate_initial_ages(parameters)
  mockery::expect_args(
    mock_exp,
    1,
    4,
    1 / parameters$average_age
  )
})

test_that('calculate_initial_ages calculates the correct proportions for a custom demographic', {
  expect_true(FALSE)
})
