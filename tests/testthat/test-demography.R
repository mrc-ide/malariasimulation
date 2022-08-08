test_that('calculate_initial_ages defaults to an exponential distribution', {
  parameters <- get_parameters(list(human_population = 4))
  mock_exp <- mockery::mock(seq(4))
  mockery::stub(calculate_initial_ages, 'rexp', mock_exp)
  ages <- calculate_initial_ages(parameters)
  mockery::expect_args(
    mock_exp,
    1,
    parameters$human_population,
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
  parameters <- get_parameters(list(human_population = 4))
  ages <- c(50, 100) * 365
  deathrates <- c(.5, .75)
  parameters <- set_demography(
    parameters,
    agegroups = ages,
    timesteps = 0,
    deathrates = matrix(deathrates, nrow = 1, ncol = 2)
  )
  mock_groups <- mockery::mock(c(2, 1, 2, 1))
  mock_rtexp <- mockery::mock(c(25 * 365, 30 * 365), c(25 * 365, 30 * 365))
  mockery::stub(calculate_initial_ages, 'sample.int', mock_groups)
  mockery::stub(calculate_initial_ages, 'rtexp', mock_rtexp)
  mockery::stub(
    calculate_initial_ages,
    'get_equilibrium_population',
    mockery::mock(c(3, 1))
  )
  ages <- calculate_initial_ages(parameters)
  mockery::expect_args(mock_groups, 1, 2, 4, replace = TRUE, prob = c(3, 1))
  mockery::expect_args(mock_rtexp, 1, 2, .5, 50 * 365)
  mockery::expect_args(mock_rtexp, 2, 2, .75, 50 * 365)
  expect_setequal(ages, c(25 * 365, 75 * 365, 30 * 365, 80 * 365))
})

test_that('match_timestep always gives 0 for constant demography', {
  expect_equal(match_timestep(c(0), 0), 1)
  expect_equal(match_timestep(c(0), 1), 1)
  expect_equal(match_timestep(c(0), 50), 1)
})

test_that('match_timestep works on the boundaries', {
  expect_equal(match_timestep(c(0, 50, 100), 0), 1)
  expect_equal(match_timestep(c(0, 50, 100), 1), 1)
  expect_equal(match_timestep(c(0, 50, 100), 49), 1)
  expect_equal(match_timestep(c(0, 50, 100), 50), 2)
  expect_equal(match_timestep(c(0, 50, 100), 99), 2)
  expect_equal(match_timestep(c(0, 50, 100), 100), 3)
  expect_equal(match_timestep(c(0, 50, 100), 101), 3)
})