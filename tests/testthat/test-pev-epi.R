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
  variables$vaccinated_timestep <- mock_integer(
    c(50, -1, -1, 4*365, -1)
  )

  events$pev_epi_doses <- lapply(events$pev_epi_doses, mock_event)

  process <- create_epi_pev_process(
    variables,
    events,
    parameters,
    get_correlation_parameters(parameters),
    parameters$pev_epi_coverages,
    parameters$pev_epi_timesteps
  )

  mockery::stub(
    process,
    'sample_intervention',
    mockery::mock(c(TRUE, TRUE, FALSE))
  )

  process(timestep)

  mockery::expect_args(
    events$pev_epi_doses[[1]]$schedule,
    1,
    c(1, 2),
    parameters$pev_doses[[1]]
  )

  mockery::expect_args(
    events$pev_epi_doses[[2]]$schedule,
    1,
    c(1, 2),
    parameters$pev_doses[[2]]
  )

  mockery::expect_args(
    events$pev_epi_doses[[3]]$schedule,
    1,
    c(1, 2),
    parameters$pev_doses[[3]]
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
