test_that('pev epi strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_pev_epi(
    parameters,
    profile = rtss_profile,
    coverages = c(0.1, 0.8),
    timesteps = c(10, 100),
    min_wait = 0,
    age = 5 * 30,
    booster_timestep = c(18, 36) * 30,
    booster_coverage = c(.9, .8),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  expect_equal(parameters$pev, TRUE)
  expect_equal(parameters$pev_epi_coverages, c(0.1, 0.8))
  expect_equal(parameters$pev_epi_timesteps, c(10, 100))
  expect_equal(parameters$pev_epi_age, 5 * 30)
  expect_equal(parameters$pev_epi_min_wait, 0)
  expect_equal(parameters$pev_epi_booster_timestep, c(18, 36) * 30)
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
      booster_timestep = c(18, 36) * 30,
      booster_coverage = c(.9, .8),
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
      booster_timestep = c(18, 36) * 30,
      booster_coverage = c(.9, .8),
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    ), "all(coverages >= 0) && all(coverages <= 1) is not TRUE",
    fixed = TRUE
  )
})

test_that('pev epi fails pre-emptively with unaligned booster parameters', {
  parameters <- get_parameters()
  expect_error(
    set_pev_epi(
      profile = rtss_profile,
      coverages = c(0.1, 0.8),
      timesteps = c(10, 100),
      min_wait = 0,
      age = 5 * 30,
      booster_timestep = c(18, 36) * 30,
      booster_coverage = .9,
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    )
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
    booster_timestep = c(18, 36) * 30,
    booster_coverage = c(.9, .8),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$birth <- individual::IntegerVariable$new(
    -c(18, 18, 2.9, 18, 18) * 365 + timestep
  )
  variables$pev_timestep <- mock_integer(
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
    variables$pev_timestep$queue_update_mock(),
    1,
    timestep,
    c(1, 2)
  )

  mockery::expect_args(
    variables$pev_profile$queue_update_mock(),
    1,
    -1,
    c(1, 2)
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
    booster_timestep = c(3, 12) * 30,
    booster_coverage = c(.9, .8),
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
    booster_timestep = c(3, 12) * 30,
    booster_coverage = c(.9, .8),
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
