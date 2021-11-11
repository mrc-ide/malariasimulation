test_that('calculate_initial_ages defaults to an exponential distribution', {
  parameters <- get_parameters(list(max_human_population = 4))
  mock_exp <- mockery::mock(seq(4))
  mockery::stub(calculate_initial_ages, 'rexp', mock_exp)
  ages <- calculate_initial_ages(parameters)
  mockery::expect_args(
    mock_exp,
    1,
    parameters$max_human_population,
    1 / parameters$average_age
  )
})

test_that('find_birthrates is consistent with get_equilibrium_population', {
  pops <- c(100, 10000, 100000)
  age_high <- c(50, 100)
  deathrates <- c(.5, .75)
  birthrates <- find_birthrates(pops, age_high, deathrates)
  actual <- vnapply(
    birthrates,
    function(birthrate) {
      sum(get_equilibrium_population(
        age_high = age_high,
        birthrate = birthrate,
        deathrates = deathrates
      ))
    }
  )
  expect_equal(actual, pops)
})


test_that('calculate_initial_ages calculates truncated exp custom demographic', {
  parameters <- get_parameters()
  parameters <- set_demography(
    parameters,
    agegroups = c(50, 100) * 365,
    timesteps = 1,
    birthrates = 50,
    deathrates = matrix(c(.5, .75), nrow=1, ncol=2)
  )
  parameters$max_human_population <- 2
  mock_groups <- mockery::mock(c(2, 1))
  mock_rtexp <- mockery::mock(25 * 365, 30 * 365)
  mockery::stub(calculate_initial_ages, 'sample.int', mock_groups)
  mockery::stub(calculate_initial_ages, 'rtexp', mock_rtexp)
  ages <- calculate_initial_ages(parameters)
  expect_setequal(ages, c(25 * 365, 80 * 365))
})
