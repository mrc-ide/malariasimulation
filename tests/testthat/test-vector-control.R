test_that('set_bednets validates coverages', {
  parameters <- get_parameters()
  expect_error(
    set_bednets(
      parameters,
      timesteps = c(5, 50),
      coverages = c(.5),
      retention = 40,
      dn0 = matrix(c(.533, .533), nrow=2, ncol=1),
      rn = matrix(c(.56, .56), nrow=2, ncol=1),
      rnm = matrix(c(.24, .24), nrow=2, ncol=1),
      gamman = c(963.6, 963.6)
    )
  )
})

test_that('set_bednets validates matrices', {
  parameters <- get_parameters()
  parameters <- set_species(parameters, list(gamb_params, fun_params), c(.1, .9))
  expect_error(
    set_bednets(
      parameters,
      timesteps = c(5, 50),
      coverages = c(.5, .9),
      retention = 40,
      dn0 = matrix(c(.533, .533), nrow=2, ncol=1),
      rn = matrix(c(.56, .56), nrow=2, ncol=1),
      rnm = matrix(c(.24, .24), nrow=2, ncol=1),
      gamman = c(963.6, 963.6)
    )
  )
})

test_that('set_bednets sets parameters', {
  parameters <- get_parameters()
  parameters <- set_bednets(
    parameters,
    timesteps = c(5, 50),
    coverages = c(.5, .9),
    retention = 40,
    dn0 = matrix(c(.533, .533), nrow=2, ncol=1),
    rn = matrix(c(.56, .56), nrow=2, ncol=1),
    rnm = matrix(c(.24, .24), nrow=2, ncol=1),
    gamman = c(963.6, 963.6)
  )
  expect_true(parameters$bednets)
  expect_equal(parameters$bednet_timesteps, c(5, 50))
  expect_equal(parameters$bednet_coverages, c(.5, .9))
  expect_equal(parameters$bednet_retention, 40)
})

test_that('set_spraying validates parameters', {
  parameters <- get_parameters()
  expect_error(
    set_spraying(
      parameters,
      timesteps = c(5, 50),
      coverages = c(.5, .9),
      ls_theta = matrix(c(2.025, 2.025), nrow=2, ncol=1),
      ls_gamma = matrix(c(-0.009, -0.009), nrow=2, ncol=1),
      ks_theta = matrix(c(-2.222, -2.222), nrow=2, ncol=1),
      ks_gamma = matrix(c(0.008, 0.008), nrow=2, ncol=1),
      ms_theta = matrix(c(-1.232, -1.232), nrow=2, ncol=1),
      ms_gamma = matrix(c(-0.009), nrow=1, ncol=1)
    )
  )
})

test_that('set_spraying sets parameters', {
  parameters <- get_parameters()
  parameters <- set_spraying(
    parameters,
    timesteps = c(5, 50),
    coverages = c(.5, .9),
    ls_theta = matrix(c(2.025, 2.025), nrow=2, ncol=1),
    ls_gamma = matrix(c(-0.009, -0.009), nrow=2, ncol=1),
    ks_theta = matrix(c(-2.222, -2.222), nrow=2, ncol=1),
    ks_gamma = matrix(c(0.008, 0.008), nrow=2, ncol=1),
    ms_theta = matrix(c(-1.232, -1.232), nrow=2, ncol=1),
    ms_gamma = matrix(c(-0.009, -0.009), nrow=2, ncol=1)
  )
  expect_true(parameters$spraying)
  expect_equal(parameters$spraying_timesteps, c(5, 50))
  expect_equal(parameters$spraying_coverages, c(.5, .9))
})

test_that('distribute_bednets process sets net_time correctly', {
  timestep <- 50
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_bednets(
    parameters,
    timesteps = c(5, 50),
    coverages = c(.5, .9),
    retention = 40,
    dn0 = matrix(c(.533, .533), nrow=2, ncol=1),
    rn = matrix(c(.56, .56), nrow=2, ncol=1),
    rnm = matrix(c(.24, .24), nrow=2, ncol=1),
    gamman = c(963.6, 963.6)
  )
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
  parameters <- set_bednets(
    parameters,
    timesteps = c(5, 50),
    coverages = c(.5, .9),
    retention = 40,
    dn0 = matrix(c(.533, .533), nrow=2, ncol=1),
    rn = matrix(c(.56, .56), nrow=2, ncol=1),
    rnm = matrix(c(.24, .24), nrow=2, ncol=1),
    gamman = c(963.6, 963.6)
  )
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
  parameters <- set_spraying(
    parameters,
    timesteps = c(5, 50),
    coverages = c(.5, .9),
    ls_theta = matrix(c(2.025, 2.025), nrow=2, ncol=1),
    ls_gamma = matrix(c(-0.009, -0.009), nrow=2, ncol=1),
    ks_theta = matrix(c(-2.222, -2.222), nrow=2, ncol=1),
    ks_gamma = matrix(c(0.008, 0.008), nrow=2, ncol=1),
    ms_theta = matrix(c(-1.232, -1.232), nrow=2, ncol=1),
    ms_gamma = matrix(c(-0.009, -0.009), nrow=2, ncol=1)
  )
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
  parameters <- get_parameters()
  parameters <- set_bednets(
    parameters,
    timesteps = c(5, 50, 100),
    coverages = c(.5, .9, .2),
    retention = 40,
    dn0 = matrix(rep(.533, 3), nrow=3, ncol=1),
    rn = matrix(rep(.56, 3), nrow=3, ncol=1),
    rnm = matrix(rep(.24, 3), nrow=3, ncol=1),
    gamman = rep(25, 3)
  )
  variables <- create_variables(parameters)
  variables$net_time <- individual::DoubleVariable$new(
    c(-1, 5, 50, 100)
  )
  variables$spray_time <- individual::DoubleVariable$new(rep(-1, 4))

  expect_equal(
    prob_bitten(timestep, variables, 1, parameters),
    list(
      prob_bitten_survives = c(1, 0.7797801, 0.6978752, 0.0709500),
      prob_bitten = c(1, 0.7797801, 0.6978752, 0.0709500),
      prob_repelled = c(0, 0.2100848, 0.2408112, 0.4760000)
    ),
    tolerance = 1e-5
  )
})

test_that('prob_bitten correctly calculates spraying only probabilities', {
  timestep <- 100
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_spraying(
    parameters,
    timesteps = c(5, 50, 100),
    coverages = c(.5, .9, .2),
    ls_theta = matrix(rep(2.025, 3), nrow=3, ncol=1),
    ls_gamma = matrix(rep(-0.009, 3), nrow=3, ncol=1),
    ks_theta = matrix(rep(-2.222, 3), nrow=3, ncol=1),
    ks_gamma = matrix(rep(0.008, 3), nrow=3, ncol=1),
    ms_theta = matrix(rep(-1.232, 3), nrow=3, ncol=1),
    ms_gamma = matrix(rep(-0.009, 3), nrow=3, ncol=1)
  )
  variables <- create_variables(parameters)

  variables$net_time <- individual::IntegerVariable$new(rep(-1, 4))
  variables$spray_time <- individual::IntegerVariable$new(
    c(-1, 5, 50, 100)
  )

  expect_equal(
    prob_bitten(timestep, variables, 1, parameters),
    list(
      prob_bitten_survives = c(1, 0.2216652, 0.1833448, 0.1506359),
      prob_bitten = c(1, 0.8268157, 0.8101383, 0.7688352),
      prob_repelled = c(0, 0.1731843, 0.1898617, 0.2311648)
    ),
    tolerance = 1e-5
  )
})

test_that('prob_bitten correctly combines spraying and net probabilities', {
  timestep <- 100
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_bednets(
    parameters,
    timesteps = c(5, 50, 100),
    coverages = c(.5, .9, .2),
    retention = 40,
    dn0 = matrix(rep(.533, 3), nrow=3, ncol=1),
    rn = matrix(rep(.56, 3), nrow=3, ncol=1),
    rnm = matrix(rep(.24, 3), nrow=3, ncol=1),
    gamman = rep(25, 3)
  )
  parameters <- set_spraying(
    parameters,
    timesteps = c(5, 50, 100),
    coverages = c(.5, .9, .2),
    ls_theta = matrix(rep(2.025, 3), nrow=3, ncol=1),
    ls_gamma = matrix(rep(-0.009, 3), nrow=3, ncol=1),
    ks_theta = matrix(rep(-2.222, 3), nrow=3, ncol=1),
    ks_gamma = matrix(rep(0.008, 3), nrow=3, ncol=1),
    ms_theta = matrix(rep(-1.232, 3), nrow=3, ncol=1),
    ms_gamma = matrix(rep(-0.009, 3), nrow=3, ncol=1)
  )
  variables <- create_variables(parameters)
  variables$net_time <- individual::IntegerVariable$new(
    c(100, 50, 5, -1)
  )
  variables$spray_time <- individual::IntegerVariable$new(
    c(-1, 5, 50, 100)
  )

  expect_equal(
    prob_bitten(timestep, variables, 1, parameters),
    list(
      prob_bitten_survives = c(0.0709500, 0.1808229, 0.1629512, 0.1506359),
      prob_bitten = c(0.0709500, 0.5828278, 0.6363754, 0.7688352),
      prob_repelled = c(0.4760000, 0.3676569, 0.3556276, 0.2311648)
    ),
    tolerance=1e-4
  )
})
