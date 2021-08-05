test_that('RTS,S strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_mass_rtss(
    parameters,
    timesteps = 10,
    coverages = 0.8,
    min_wait = 0,
    min_ages = 5 * 30,
    max_ages = 17 * 30,
    boosters = c(18, 36) * 30,
    booster_coverage = c(.9, .8)
  )
  expect_equal(parameters$rtss, TRUE)
  expect_equal(parameters$rtss_mass_timesteps, 10)
  expect_equal(parameters$rtss_mass_coverages, .8)
  expect_equal(parameters$rtss_mass_min_ages, 5 * 30)
  expect_equal(parameters$rtss_mass_max_ages, 17 * 30)
  expect_equal(parameters$rtss_mass_boosters, c(18, 36) * 30)
})

test_that('RTS,S fails pre-emptively', {
  parameters <- get_parameters()
  expect_error(
    set_mass_rtss(
      parameters,
      timesteps = 10,
      coverages = 0.8,
      min_wait = 0,
      min_ages = 5 * 30,
      max_ages = 17 * 30,
      boosters = c(18, 36) * 30,
      coverage = 0.8,
      booster_coverage = .9
    )
  )
  expect_error(
    set_rtss(
      parameters,
      timesteps = 10,
      coverages = 0.8,
      min_ages = c(0, 5 * 30),
      max_ages = 17 * 30,
      boosters = c(18, 36) * 30,
      booster_coverage = .9
    )
  )
})

test_that('Infection considers vaccine efficacy', {
  timestep <- 100
  parameters <- get_parameters()
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  variables$state <- individual::CategoricalVariable$new(
    c('S', 'A', 'D', 'U', 'Tr'),
    rep('S', 4)
  )
  variables$birth <- individual::IntegerVariable$new(
    -c(8, 2.9, 3.2, 18.4) * 365 - 100
  )
  variables$rtss_vaccinated <- individual::IntegerVariable$new(
    c(-1, -1, 50, 50)
  )
  variables$rtss_boosted <- individual::IntegerVariable$new(
    c(-1, -1, -1, 50 + 30)
  )
  variables$rtss_cs <- individual::DoubleVariable$new(exp(
    c(rep(parameters$rtss_cs[[1]], 3), parameters$rtss_cs_boost[[1]])
  ))
  variables$rtss_rho <- individual::DoubleVariable$new(invlogit(
    c(rep(parameters$rtss_rho[[1]], 3), parameters$rtss_rho_boost[[1]])
  ))
  variables$rtss_ds <- individual::DoubleVariable$new(
    exp(rep(parameters$rtss_ds[[1]], 4))
  )
  variables$rtss_dl <- individual::DoubleVariable$new(
    exp(rep(parameters$rtss_dl[[1]], 4))
  )
  variables$drug <- individual::IntegerVariable$new(
    c(-1, -1, -1, -1)
  )
  variables$ib <- individual::IntegerVariable$new(
    rep(.2, 4)
  )

  bernoulli_mock <- mockery::mock(c(1, 2))
  mockery::stub(calculate_infections, 'bernoulli_multi_p', bernoulli_mock)
  calculate_infections(
    variables,
    bitten_humans = individual::Bitset$new(4)$insert(seq(4)),
    parameters,
    timestep
  )

  expect_equal(
    mockery::mock_args(bernoulli_mock)[[1]][[1]],
    c(0.590, 0.590, 0.215, 0.244),
    tolerance=1e-3
  )
})

test_that('RTS,S vaccinations update vaccination time and schedule boosters', {
  timestep <- 100
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_mass_rtss(
    parameters,
    timesteps = c(100, 100 + 365),
    coverages = rep(0.8, 2),
    min_wait = 0,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    boosters = c(18, 36) * 30,
    booster_coverage = c(.9, .8)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$birth <- individual::IntegerVariable$new(
    -c(18.3, 8, 2.9, 3.2, 18.4) * 365 + 100
  )
  variables$rtss_vaccinated <- mock_integer(
    c(-1, -1, -1, 50, 50)
  )

  events$rtss_mass_vaccination <- mock_event(events$rtss_mass_vaccination)
  events$rtss_mass_doses <- lapply(events$rtss_mass_doses, mock_event)

  renderer <- individual::Render$new(100)

  listener <- create_rtss_mass_listener(
    variables,
    events,
    parameters,
    get_correlation_parameters(parameters)
  )

  mockery::stub(
    listener,
    'sample_intervention',
    mockery::mock(c(TRUE, TRUE, FALSE, FALSE))
  )

  listener(timestep)

  mockery::expect_args(
    events$rtss_mass_doses[[1]]$schedule,
    1,
    c(1, 3),
    0
  )

  mockery::expect_args(
    events$rtss_mass_doses[[2]]$schedule,
    1,
    c(1, 3),
    7
  )

  mockery::expect_args(
    events$rtss_mass_doses[[3]]$schedule,
    1,
    c(1, 3),
    14
  )

  mockery::expect_args(
    events$rtss_mass_vaccination$schedule,
    1,
    365
  )
})

test_that('RTS,S boosters update antibody params and reschedule correctly', {
  parameters <- get_parameters()
  parameters <- set_mass_rtss(
    parameters,
    timesteps = c(50, 50 + 365),
    coverages = rep(0.8, 2),
    min_wait = 0,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    boosters = c(1, 6) * 30,
    booster_coverage = c(1, 1)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  timestep <- 50 + 30

  variables$birth <- individual::IntegerVariable$new(
    -c(2.9, 3.2, 18.4) * 365 + 100
  )
  variables$rtss_vaccinated <- mock_double(
    c(50, 50, 50)
  )
  variables$rtss_boosted <- mock_double(
    c(-1, -1, -1)
  )
  variables$rtss_cs <- mock_double(
    c(0, 0, 0)
  )
  variables$rtss_rho <- mock_double(
    c(0, 0, 0)
  )
  events$rtss_booster <- mock_event(events$rtss_booster)

  listener <- create_rtss_booster_listener(
    variables,
    events,
    parameters,
    parameters$rtss_mass_boosters,
    parameters$rtss_mass_booster_coverage
  )

  mockery::stub(listener, 'bernoulli', mockery::mock(c(1, 2, 3)))
  mockery::stub(listener, 'rnorm', mockery::mock(c(0, 0, 0), c(0, 0, 0)))
  listener(timestep, individual::Bitset$new(3)$insert(c(1, 2, 3)))

  expect_bitset_update(
    variables$rtss_cs$queue_update_mock(),
    rep(exp(parameters$rtss_cs_boost[[1]]), 3),
    c(1, 2, 3)
  )

  expect_bitset_update(
    variables$rtss_rho$queue_update_mock(),
    rep(invlogit(parameters$rtss_rho_boost[[1]]), 3),
    c(1, 2, 3)
  )

  expect_bitset_update(
    variables$rtss_boosted$queue_update_mock(),
    timestep,
    c(1, 2, 3)
  )

  mockery::expect_args(
    events$rtss_booster$schedule,
    1,
    c(1, 2, 3),
    5 * 30
  )
})

test_that('RTS,S booster coverages sample subpopulations correctly', {
  parameters <- get_parameters()
  parameters <- set_mass_rtss(
    parameters,
    timesteps = 50,
    coverages = 0.8,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    min_wait = 0,
    boosters = c(1, 6) * 30,
    booster_coverage = c(.9, .8)
  )

  timestep <- 50 + 30

  events <- create_events(parameters)
  events$rtss_booster <- mock_event(events$rtss_booster)
  variables <- create_variables(parameters)

  variables$birth <- individual::IntegerVariable$new(
    -c(2.9, 3.2, 18.4) * 365 + 100
  )
  variables$rtss_vaccinated <- mock_double(
    c(50, 50, 50)
  )
  variables$rtss_boosted <- mock_double(
    c(-1, -1, -1)
  )
  variables$rtss_cs <- mock_double(
    c(0, 0, 0)
  )
  variables$rtss_rho <- mock_double(
    c(0, 0, 0)
  )

  listener <- create_rtss_booster_listener(
    variables,
    events,
    parameters,
    parameters$rtss_mass_boosters,
    parameters$rtss_mass_booster_coverage
  )

  bernoulli_mock = mockery::mock(c(2, 3))

  mockery::stub(listener, 'bernoulli', bernoulli_mock)
  mockery::stub(listener, 'rnorm', mockery::mock(c(0, 0, 0), c(0, 0, 0)))
  listener(timestep, individual::Bitset$new(3)$insert(c(1, 2, 3)))

  mockery::expect_args(
    bernoulli_mock,
    1,
    3, #number of humans
    .8 #booster coverage for the second booster
  )

  expect_bitset_update(
    variables$rtss_cs$queue_update_mock(),
    rep(exp(parameters$rtss_cs_boost[[1]]), 3),
    c(1, 2, 3)
  )

  expect_bitset_update(
    variables$rtss_rho$queue_update_mock(),
    rep(invlogit(parameters$rtss_rho_boost[[1]]), 3),
    c(1, 2, 3)
  )

  expect_bitset_update(
    variables$rtss_boosted$queue_update_mock(),
    timestep,
    c(1, 2, 3)
  )

  mockery::expect_args(
    events$rtss_booster$schedule,
    1,
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
