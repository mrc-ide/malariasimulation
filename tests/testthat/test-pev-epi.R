test_that('pev epi strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_pev_epi(
    parameters,
    profile = rtss_profile,
    coverages = c(0.1, 0.8),
    timesteps = c(10, 100),
    min_wait = 0,
    age = 5 * 30,
    booster_spacing = c(18, 36) * 30,
    booster_coverage = matrix(c(.9, .8, .9, .8), nrow=2, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  expect_equal(parameters$pev, TRUE)
  expect_equal(parameters$pev_epi_coverages, c(0.1, 0.8))
  expect_equal(parameters$pev_epi_timesteps, c(10, 100))
  expect_equal(parameters$pev_epi_age, 5 * 30)
  expect_equal(parameters$pev_epi_min_wait, 0)
  expect_equal(parameters$pev_epi_booster_spacing, c(18, 36) * 30)
  expect_equal(parameters$pev_profiles, list(rtss_profile, rtss_booster_profile, rtss_booster_profile))
  expect_equal(parameters$pev_epi_profile_indices, seq(3))

  expect_error(
    parameters <- set_pev_epi(
      parameters,
      profile = rtss_profile,
      coverages = -1, # less than 0
      timesteps = 10,
      min_wait = 0,
      age = 5 * 30,
      booster_spacing = c(18, 36) * 30,
      booster_coverage = matrix(c(.9, .8, .9, .8), nrow=2, ncol=2),
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    ), "all(coverages >= 0) && all(coverages <= 1) is not TRUE",
    fixed = TRUE
  )

  expect_error(
    parameters <- set_pev_epi(
      parameters,
      profile = rtss_profile,
      coverages = 1.5, # greater than 1
      timesteps = 10,
      min_wait = 0,
      age = 5 * 30,
      booster_spacing = c(18, 36) * 30,
      booster_coverage = matrix(c(.9, .8, .9, .8), nrow=2, ncol=2),
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    ), "all(coverages >= 0) && all(coverages <= 1) is not TRUE",
    fixed = TRUE
  )
})

test_that('set_pev_epi checks booster coverage matrix shape', {
  parameters <- get_parameters()
  expect_error(
    parameters <- set_pev_epi(
      parameters,
      profile = rtss_profile,
      coverages = c(0.1, 0.8),
      timesteps = c(10, 100),
      min_wait = 0,
      age = 5 * 30,
      booster_spacing = c(18, 36) * 30,
      booster_coverage = matrix(c(.9, .8), nrow=2, ncol=1),
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    ),
    'booster_spacing, booster_coverage and booster_profile do not align',
    fixed = TRUE
  )
})

test_that('pev epi targets correct age and respects min_wait', {
  timestep <- 5*365 
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_pev_epi(
    parameters,
    profile = rtss_profile,
    timesteps = 10,
    coverages = 0.8,
    min_wait = 2*365,
    age = 18 * 365,
    booster_spacing = c(18, 36) * 30,
    booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$birth <- individual::IntegerVariable$new(
    -c(18, 18, 2.9, 18, 18) * 365 + timestep
  )
  variables$last_pev_timestep <- mock_integer(
    c(50, -1, -1, 4*365, -1)
  )
  variables$pev_profile <- mock_integer(
    c(1, -1, -1, 1, -1)
  )

  correlations <- get_correlation_parameters(parameters)
  process <- create_epi_pev_process(
    variables,
    events,
    parameters,
    correlations,
    parameters$pev_epi_coverages,
    parameters$pev_epi_timesteps
  )

  sample_mock <- mockery::mock(c(TRUE, TRUE, FALSE))
  mockery::stub(
    process,
    'sample_intervention',
    sample_mock
  )

  process(timestep)

  mockery::expect_args(
    sample_mock,
    1,
    c(1, 2, 5),
    'pev',
    .8,
    correlations
  )

  mockery::expect_args(
    variables$last_pev_timestep$queue_update_mock(),
    1,
    timestep,
    c(1, 2)
  )

})

test_that('EPI ignores individuals scheduled for mass vaccination', {
  timestep <- 100
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = c(50, 100),
    coverages = rep(0.8, 2),
    min_wait = 0,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    booster_spacing = c(18, 36) * 30,
    booster_coverage = matrix(c(.9, .8, .9, .8), nrow=2, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  parameters <- set_pev_epi(
    parameters,
    profile = rtss_profile,
    timesteps = 10,
    coverages = 0.8,
    min_wait = 0,
    age = 18 * 365,
    booster_spacing = c(18, 36) * 30,
    booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$birth <- individual::IntegerVariable$new(
    -c(18, 8, 2.9, 3.2, 18) * 365 + 100
  )
  variables$pev_timestep <- mock_integer(
    c(-1, -1, -1, 50, 50)
  )

  correlations <- get_correlation_parameters(parameters)

  process <- create_epi_pev_process(
    variables,
    events,
    parameters,
    correlations,
    parameters$pev_epi_coverages,
    parameters$pev_epi_timesteps
  )

  sample_mock <- mockery::mock(c(TRUE))

  mockery::stub(
    process,
    'sample_intervention',
    sample_mock
  )
  
  # schedule id #1 for epi vaccination
  events$mass_pev_doses[[1]]$schedule(1, 0)

  process(timestep)

  mockery::expect_args(
    sample_mock,
    1,
    5,
    'pev',
    .8,
    correlations
  )
})


test_that('pev EPI respects min_wait when scheduling seasonal boosters', {
  timestep <- 5 * 365 
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_pev_epi(
    parameters,
    profile = rtss_profile,
    timesteps = 10,
    coverages = 0.8,
    min_wait = 6 * 30,
    age = 18 * 365,
    booster_spacing = c(3, 12) * 30,
    booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile),
    seasonal_boosters = TRUE
  )
  events <- create_events(parameters)

  booster_event <- mock_event(events$pev_epi_boosters[[1]])

  listener <- create_seasonal_booster_scheduler(
    booster_event = booster_event,
    booster_delay = 3 * 30,
    parameters
  )

  target <- individual::Bitset$new(5)$insert(seq(5))

  listener(timestep, target)

  expect_bitset_schedule(
    booster_event$schedule,
    seq(5),
    365 + 3 * 30
  )
})

test_that('pev EPI schedules for the following year with seasonal boosters', {
  timestep <- 5 * 365 + 4 * 30
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_pev_epi(
    parameters,
    profile = rtss_profile,
    timesteps = 10,
    coverages = 0.8,
    min_wait = 6 * 30,
    age = 18 * 365,
    booster_spacing = c(3, 12) * 30,
    booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile),
    seasonal_boosters = TRUE
  )
  events <- create_events(parameters)

  booster_event <- mock_event(events$pev_epi_boosters[[1]])

  listener <- create_seasonal_booster_scheduler(
    booster_event = booster_event,
    booster_delay = 3 * 30,
    parameters
  )

  target <- individual::Bitset$new(5)$insert(seq(5))

  listener(timestep, target)

  expect_bitset_schedule(
    booster_event$schedule,
    seq(5),
    335
  )
})
