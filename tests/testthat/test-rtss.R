test_that('RTS,S strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_rtss(
    parameters,
    start = 10,
    end = 100,
    frequency = 365,
    min_ages = 5 * 30,
    max_ages = 17 * 30,
    boosters = c(18, 36) * 30,
    coverage = 0.8,
    booster_coverage = c(.9, .8)
  )
  expect_equal(parameters$rtss, TRUE)
  expect_equal(parameters$rtss_start, 10)
  expect_equal(parameters$rtss_end, 100)
  expect_equal(parameters$rtss_min_ages, 5 * 30)
  expect_equal(parameters$rtss_max_ages, 17 * 30)
  expect_equal(parameters$rtss_boosters, c(18, 36) * 30)
  expect_equal(parameters$rtss_coverage, .8)
})

test_that('RTS,S fails pre-emptively', {
  parameters <- get_parameters()
  expect_error(
    set_rtss(
      parameters,
      start = 100,
      end = 100,
      frequency = 365,
      min_ages = 5 * 30,
      max_ages = 17 * 30,
      boosters = c(18, 36) * 30,
      coverage = 0.8,
      booster_coverage = c(.9)
    ),
    '*'
  )
  expect_error(
    set_rtss(
      parameters,
      start = 100,
      end = 100,
      frequency = 365,
      min_ages = c(0, 5 * 30),
      max_ages = 17 * 30,
      boosters = c(18, 36) * 30,
      coverage = 0.8,
      booster_coverage = c(.9)
    ),
    '*'
  )
  expect_error(
    set_rtss(
      parameters,
      start = 10,
      end = 100,
      frequency = 365,
      min_ages = 5 * 30,
      max_ages = 17 * 30,
      boosters = c(18, 36) * 30,
      coverage = 0.8,
      booster_coverage = c(.9)
    ),
    '*'
  )
})

test_that('Infection considers vaccine efficacy', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        S = seq(4),
        birth = -c(8, 2.9, 3.2, 18.4) * 365 - 100,
        rtss_vaccinated = c(-1, -1, 50, 50),
        rtss_boosted = c(-1, -1, -1, 50 + 30),
        rtss_cs = exp(
          c(rep(parameters$rtss_cs[[1]], 3), parameters$rtss_cs_boost[[1]])
        ),
        rtss_rho = invlogit(
          c(rep(parameters$rtss_rho[[1]], 3), parameters$rtss_rho_boost[[1]])
        ),
        rtss_ds = exp(rep(parameters$rtss_ds[[1]], 4)),
        rtss_dl = exp(rep(parameters$rtss_dl[[1]], 4)),
        drug = c(-1, -1, -1, -1)
      )
    ),
    parameters = parameters,
    timestep = 100
  )

  bernoulli_mock <- mockery::mock(c(TRUE, TRUE, FALSE, FALSE))
  with_mock(
    'malariasimulation:::bernoulli_multi_p' = bernoulli_mock,
    calculate_infections(
      api,
      individuals$human,
      states,
      variables,
      seq(4),
      rep(.2, 4)
    )
  )

  expect_equal(
    mockery::mock_args(bernoulli_mock)[[1]][[2]],
    c(0.590, 0.590, 0.215, 0.244),
    tolerance=1e-3
  )
})

test_that('RTS,S vaccinations update vaccination time and schedule boosters', {
  parameters <- get_parameters()
  parameters <- set_rtss(
    parameters,
    start = 50,
    end = 100 + 365,
    frequency = 365,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    boosters = c(18, 36) * 30,
    coverage = 0.8,
    booster_coverage = c(.9, .8)
  )
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        birth = -c(18.3, 8, 2.9, 3.2, 18.4) * 365 + 100,
        rtss_vaccinated = c(-1, -1, -1, 50, 50)
      )
    ),
    parameters = parameters,
    timestep = 100
  )

  listener <- create_rtss_vaccination_listener(
    individuals$human,
    variables,
    events,
    parameters,
    get_correlation_parameters(parameters)
  )

  bernoulli_mock = mockery::mock(2)

  with_mock(
    'malariasimulation:::bernoulli' = bernoulli_mock,
    'malariasimulation:::sample_intervention' = mockery::mock(TRUE, TRUE, FALSE),
    listener(api, c(1))
  )

  mockery::expect_args(
    bernoulli_mock,
    1, # second call
    2, # n vaccinated
    .9 # first booster coverage
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$rtss_vaccinated,
    100,
    c(1, 3)
  )

  mockery::expect_args(
    api$schedule,
    1,
    events$rtss_vaccination,
    c(1),
    365
  )
  mockery::expect_args(
    api$schedule,
    2,
    events$rtss_booster,
    3,
    18 * 30
  )
})

test_that('RTS,S boosters update antibody params and reschedule correctly', {
  parameters <- get_parameters()
  parameters <- set_rtss(
    parameters,
    start = 50,
    end = 200,
    frequency = 365,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    boosters = c(1, 6) * 30,
    coverage = 0.8,
    booster_coverage = c(1, 1)
  )
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        birth = -c(2.9, 3.2, 18.4) * 365 + 100,
        rtss_vaccinated = c(50, 50, 50),
        rtss_boosted = c(-1, -1, -1)
      )
    ),
    parameters = parameters,
    timestep = 50 + 30
  )

  listener <- create_rtss_booster_listener(
    individuals$human,
    variables,
    events,
    parameters
  )

  with_mock(
    "malariasimulation::bernoulli" = mockery::mock(c(1, 2, 3)),
    rnorm = mockery::mock(c(0, 0, 0), c(0, 0, 0)),
    listener(api, c(1, 2, 3))
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$rtss_cs,
    rep(exp(parameters$rtss_cs_boost[[1]]), 3),
    c(1, 2, 3)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$rtss_rho,
    rep(invlogit(parameters$rtss_rho_boost[[1]]), 3),
    c(1, 2, 3)
  )

  mockery::expect_args(
    api$queue_variable_update,
    3,
    individuals$human,
    variables$rtss_boosted,
    50 + 30,
    c(1, 2, 3)
  )

  mockery::expect_args(
    api$schedule,
    1,
    events$rtss_booster,
    c(1, 2, 3),
    5 * 30
  )
})

test_that('RTS,S booster coverages sample subpopulations correctly', {
  parameters <- get_parameters()
  parameters <- set_rtss(
    parameters,
    start = 50,
    end = 200,
    frequency = 365,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    boosters = c(1, 6) * 30,
    coverage = 0.8,
    booster_coverage = c(.9, .8)
  )
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        birth = -c(2.9, 3.2, 18.4) * 365 + 100,
        rtss_vaccinated = c(50, 50, 50),
        rtss_boosted = c(-1, -1, -1)
      )
    ),
    parameters = parameters,
    timestep = 50 + 30
  )

  listener <- create_rtss_booster_listener(
    individuals$human,
    variables,
    events,
    parameters
  )

  bernoulli_mock = mockery::mock(c(2, 3))

  with_mock(
    rnorm = mockery::mock(c(0, 0, 0), c(0, 0, 0)),
    'malariasimulation:::bernoulli' = bernoulli_mock,
    listener(api, c(1, 2, 3))
  )

  mockery::expect_args(
    bernoulli_mock,
    1,
    3, #number of humans
    .8 #booster coverage for the second booster
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$rtss_cs,
    rep(exp(parameters$rtss_cs_boost[[1]]), 3),
    c(1, 2, 3)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$rtss_rho,
    rep(invlogit(parameters$rtss_rho_boost[[1]]), 3),
    c(1, 2, 3)
  )

  mockery::expect_args(
    api$queue_variable_update,
    3,
    individuals$human,
    variables$rtss_boosted,
    50 + 30,
    c(1, 2, 3)
  )

  mockery::expect_args(
    api$schedule,
    1,
    events$rtss_booster,
    c(2, 3),
    5 * 30
  )
})

test_that('RTS,S Antibodies are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_rtss_antibodies(
      c(0, 0, 10, 30),
      exp(c(6, 6, 5, 5)),
      invlogit(c(2, 2, 1, 1)),
      exp(c(3, 3, 3, 3)),
      exp(c(6, 6, 6, 6)),
      parameters
    ),
    c(403.4, 403.4, 116.1, 76.4),
    tolerance = 1e-3
  )
})

test_that('Efficacies are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_rtss_efficacy(c(400, 400, 100, 50), parameters),
    c(0.685, 0.685, 0.466, 0.349),
    tolerance = 1e-3
  )
})
