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
  timestep <- 50
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_bednets(parameters, c(5, 50), c(.5, .9), 40)
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$net_time <- mock_double(rep(0, 4))
  events$throw_away_net <- mock_event(events$throw_away_net)
  correlations <- get_correlation_parameters(parameters)
  process <- distribute_nets(
    variables,
    events$throw_away_net,
    parameters,
    correlations
  )

  target_mock <- mockery::mock(c(FALSE, FALSE, TRUE, TRUE))
  mockery::stub(process, 'sample_intervention', target_mock)
  mockery::stub(process, 'log_uniform', mockery::mock(c(3, 4)))

  process(timestep)

  mockery::expect_args(target_mock, 1, seq(4), 'bednets', .9, correlations)
  mockery::expect_args(
    variables$net_time$queue_update,
    1,
    50,
    c(3, 4)
  )
  mockery::expect_args(
    events$throw_away_net$schedule,
    1,
    c(3, 4),
    c(3, 4)
  )
})

test_that('throw_away_bednets process resets net_time correctly', {
  timestep <- 1
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_bednets(parameters, c(5, 50), c(.5, .9), 40)
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$net_time <- mock_double(rep(0, 4))
  listener <- throw_away_nets(variables)

  listener(timestep, individual::Bitset$new(4)$insert(c(2, 3)))

  expect_bitset_update(
    variables$net_time$queue_update,
    -1,
    c(2, 3)
  )
})

test_that('indoor_spraying process sets spray_time correctly', {
  timestep <- 50
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_spraying(parameters, c(5, 50), c(.5, .9))
  spray_time <- mock_double(rep(0, 4))
  correlations <- get_correlation_parameters(parameters)
  process <- indoor_spraying(
    spray_time,
    parameters,
    correlations
  )

  target_mock <- mockery::mock(c(FALSE, FALSE, TRUE, TRUE))
  mockery::stub(process, 'sample_intervention', target_mock)

  process(timestep)

  mockery::expect_args(target_mock, 1, seq(4), 'spraying', .9, correlations)
  mockery::expect_args(
    spray_time$queue_update,
    1,
    50,
    c(3, 4)
  )
})

test_that('prob_bitten defaults to 1 with no protection', {
  timestep <- 100
  parameters <- get_parameters(list(human_population = 4))
  variables <- create_variables(parameters)
  variables$net_time <- individual::DoubleVariable$new(rep(-1, 4))
  variables$spray_time <- individual::DoubleVariable$new(rep(-1, 4))

  expect_equal(
    prob_bitten(timestep, variables, 1, parameters),
    list(
      prob_bitten_survives = rep(1, 4),
      prob_bitten = rep(1, 4),
      prob_repelled = rep(0, 4)
    )
  )
})

test_that('prob_bitten correctly calculates net only probabilities', {
  timestep <- 100
  parameters <- get_parameters(list(
    bednets = TRUE,
    gamman = 25
  ))
  variables <- create_variables(parameters)
  variables$net_time <- individual::DoubleVariable$new(
    c(-1, 5, 50, 100)
  )
  variables$spray_time <- individual::DoubleVariable$new(rep(-1, 4))

  expect_equal(
    prob_bitten(timestep, variables, 1, parameters),
    list(
      prob_bitten_survives = c(1, 0.7694168, 0.6836575, 0.0272300),
      prob_bitten = c(1, 0.7694168, 0.6836575, 0.0272300),
      prob_repelled = c(0, 0.2199712, 0.2521435, 0.4984000)
    ),
    tolerance = 1e-5
  )
})

test_that('prob_bitten correctly calculates spraying only probabilities', {
  timestep <- 100
  parameters <- get_parameters(list(spraying = TRUE, gammas = 25))
  variables <- create_variables(parameters)

  variables$net_time <- individual::DoubleVariable$new(rep(-1, 4))
  variables$spray_time <- individual::DoubleVariable$new(
    c(-1, 5, 50, 100)
  )

  expect_equal(
    prob_bitten(timestep, variables, 1, parameters),
    list(
      prob_bitten_survives = c(1, 0.9780972, 0.8699070, 0.1751120),
      prob_bitten = c(1, 0.9956601, 0.9737450, 0.8060000),
      prob_repelled = c(0, 0.00433993, 0.02625504, 0.19400000)
    ),
    tolerance = 1e-5
  )
})

test_that('prob_bitten correctly combines spraying and net probabilities', {
  timestep <- 100
  parameters <- get_parameters(list(bednets = TRUE, spraying = TRUE, gamman = 25))
  variables <- create_variables(parameters)
  variables$net_time <- individual::DoubleVariable$new(
    c(100, 50, 5, -1)
  )
  variables$spray_time <- individual::DoubleVariable$new(
    c(-1, 5, 50, 100)
  )

  expect_equal(
    prob_bitten(timestep, variables, 1, parameters),
    list(
      prob_bitten_survives = c(0.0272300, 0.4631212, 0.3765613, 0.1751120),
      prob_bitten = c(0.0272300, 0.6375005, 0.6839200, 0.8060000),
      prob_repelled = c(0.4984000, 0.3028339, 0.3066950, 0.1940000)
    ),
    tolerance=1e-4
  )
})
