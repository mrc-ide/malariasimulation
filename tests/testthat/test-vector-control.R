test_that('set_bednets validates parameters', {
  parameters <- get_parameters()
  expect_error(
    set_bednets(parameters, c(5, 50), .8, 40),
    '*'
  )
})

test_that('set_bednets sets parameters', {
  parameters <- get_parameters()
  parameters <- set_bednets(parameters, c(5, 50), c(.5, .9), 40)
  expect_true(parameters$bednets)
  expect_equal(parameters$bednet_timesteps, c(5, 50))
  expect_equal(parameters$bednet_coverages, c(.5, .9))
  expect_equal(parameters$bednet_retention, 40)
})

test_that('set_spraying validates parameters', {
  parameters <- get_parameters()
  expect_error(
    set_spraying(parameters, c(5, 50), .8),
    '*'
  )
})

test_that('set_spraying sets parameters', {
  parameters <- get_parameters()
  parameters <- set_spraying(parameters, c(5, 50), c(.5, .9))
  expect_true(parameters$spraying)
  expect_equal(parameters$spraying_timesteps, c(5, 50))
  expect_equal(parameters$spraying_coverages, c(.5, .9))
})

test_that('distribute_bednets process sets net_time correctly', {
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_bednets(parameters, c(5, 50), c(.5, .9), 40)
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  process <- distribute_nets(individuals$human, variables, parameters)

  api <- mock_api(
    list(),
    timestep = 50,
    parameters = parameters
  )

  bernoulli_mock <- mockery::mock(c(3, 4))
  mockery::stub(process, 'bernoulli', bernoulli_mock)
  mockery::stub(process, 'runif', mockery::mock(c(.99999, .5)))

  process(api)

  mockery::expect_args(bernoulli_mock, 1, 4, .9)
  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$net_time,
    50,
    c(3, 4)
  )
  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$net_end_time,
    c(51, 78),
    c(3, 4)
  )
})

test_that('throw_away_bednets process resets net_time correctly', {
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_bednets(parameters, c(5, 50), c(.5, .9), 40)
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  process <- throw_away_nets(individuals$human, variables)

  api <- mock_api(
    list(
      human = list(
        net_time = c(-1, 50, 50, 50),
        net_end_time = c(-1, 55, 55, 56)
      )
    ),
    timestep = 55,
    parameters = parameters
  )

  process(api)

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$net_time,
    -1,
    c(2, 3)
  )
  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$net_end_time,
    -1,
    c(2, 3)
  )
})

test_that('indoor_spraying process sets spray_time correctly', {
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_spraying(parameters, c(5, 50), c(.5, .9))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  process <- indoor_spraying(individuals$human, variables$spray_time, parameters)

  api <- mock_api(
    list(),
    timestep = 50,
    parameters = parameters
  )

  bernoulli_mock <- mockery::mock(c(3, 4))
  mockery::stub(process, 'bernoulli', bernoulli_mock)

  process(api)

  mockery::expect_args(bernoulli_mock, 1, 4, .9)
  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$spray_time,
    50,
    c(3, 4)
  )
})

test_that('prob_bitten defaults to 1 with no protection', {
  parameters <- get_parameters(list(human_population = 4))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        net_time = c(-1, -1, -1, -1),
        spray_time = c(-1, -1, -1, -1)
      )
    ),
    timestep = 100,
    parameters = parameters
  )

  expect_equal(
    prob_bitten(individuals, variables, 1, api, parameters),
    list(
      prob_bitten_survives = rep(1, 4),
      prob_bitten = rep(1, 4),
      prob_repelled = rep(0, 4)
    )
  )
})

test_that('prob_bitten correctly calculates net only probabilities', {
  parameters <- get_parameters(list(
    bednets = TRUE,
    gamman = 25
  ))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        net_time = c(-1, 5, 50, 100),
        spray_time = c(-1, -1, -1, -1)
      )
    ),
    timestep = 100,
    parameters = parameters
  )
  expect_equal(
    prob_bitten(individuals, variables, 1, api, parameters),
    list(
      prob_bitten_survives = c(1, 0.7694168, 0.6836575, 0.0272300),
      prob_bitten = c(1, 0.7694168, 0.6836575, 0.0272300),
      prob_repelled = c(0, 0.2199712, 0.2521435, 0.4984000)
    ),
    tolerance = 1e-5
  )
})

test_that('prob_bitten correctly calculates spraying only probabilities', {
  parameters <- get_parameters(list(spraying = TRUE, gammas = 25))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        net_time = c(-1, -1, -1, -1),
        spray_time = c(-1, 5, 50, 100)
      )
    ),
    timestep = 100,
    parameters = parameters
  )

  expect_equal(
    prob_bitten(individuals, variables, 1, api, parameters),
    list(
      prob_bitten_survives = c(1, 0.9780972, 0.8699070, 0.1751120),
      prob_bitten = c(1, 0.9956601, 0.9737450, 0.8060000),
      prob_repelled = c(0, 0.00433993, 0.02625504, 0.19400000)
    ),
    tolerance = 1e-5
  )
})

test_that('prob_bitten correctly combines spraying and net probabilities', {
  parameters <- get_parameters(list(bednets = TRUE, spraying = TRUE, gamman = 25))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        net_time = c(100, 50, 5, -1),
        spray_time = c(-1, 5, 50, 100)
      )
    ),
    timestep = 100,
    parameters = parameters
  )

  expect_equal(
    prob_bitten(individuals, variables, 1, api, parameters),
    list(
      prob_bitten_survives = c(0.0272300, 0.4631212, 0.3765613, 0.1751120),
      prob_bitten = c(0.0272300, 0.6375005, 0.6839200, 0.8060000),
      prob_repelled = c(0.4984000, 0.3028339, 0.3066950, 0.1940000)
    ),
    tolerance=1e-4
  )
})
