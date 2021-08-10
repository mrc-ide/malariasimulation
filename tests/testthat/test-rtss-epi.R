test_that('RTS,S epi strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_rtss_epi(
    parameters,
    start = 10,
    end = 100,
    coverage = 0.8,
    min_wait = 0,
    age = 5 * 30,
    boosters = c(18, 36) * 30,
    booster_coverage = c(.9, .8)
  )
  expect_equal(parameters$rtss, TRUE)
  expect_equal(parameters$rtss_epi_start, 10)
  expect_equal(parameters$rtss_epi_end, 100)
  expect_equal(parameters$rtss_epi_coverage, .8)
  expect_equal(parameters$rtss_epi_age, 5 * 30)
  expect_equal(parameters$rtss_epi_min_wait, 0)
  expect_equal(parameters$rtss_epi_boosters, c(18, 36) * 30)
})

test_that('RTS,S epi fails pre-emptively', {
  parameters <- get_parameters()
  expect_error(
    set_rtss_epi(
      parameters,
      start = 10,
      end = 100,
      coverages = 0.8,
      min_wait = 0,
      min_ages = 5 * 30,
      max_ages = 17 * 30,
      boosters = c(18, 36) * 30,
      coverage = 0.8,
      booster_coverage = .9
    )
  )
})

test_that('RTS,S epi targets correct age and respects min_wait', {
  timestep <- 5*365 
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_rtss_epi(
    parameters,
    start = 10,
    end = timestep,
    coverage = 0.8,
    min_wait = 2*365,
    age = 18 * 365,
    boosters = c(18, 36) * 30,
    booster_coverage = c(.9, .8)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$birth <- individual::IntegerVariable$new(
    -c(18, 18, 2.9, 18, 18) * 365 + timestep
  )
  variables$rtss_vaccinated <- mock_integer(
    c(50, -1, -1, 4*365, -1)
  )

  events$rtss_epi_doses <- lapply(events$rtss_epi_doses, mock_event)

  process <- create_rtss_epi_process(
    variables,
    events,
    parameters,
    get_correlation_parameters(parameters)
  )

  mockery::stub(
    process,
    'sample_intervention',
    mockery::mock(c(TRUE, TRUE, FALSE))
  )

  process(timestep)

  mockery::expect_args(
    events$rtss_epi_doses[[1]]$schedule,
    1,
    c(1, 2),
    parameters$rtss_doses[[1]]
  )

  mockery::expect_args(
    events$rtss_epi_doses[[2]]$schedule,
    1,
    c(1, 2),
    parameters$rtss_doses[[2]]
  )

  mockery::expect_args(
    events$rtss_epi_doses[[3]]$schedule,
    1,
    c(1, 2),
    parameters$rtss_doses[[3]]
  )
})

test_that('RTS,S EPI respects min_wait when scheduling seasonal boosters', {
  timestep <- 5 * 365 
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_rtss_epi(
    parameters,
    start = 10,
    end = timestep,
    coverage = 0.8,
    min_wait = 6 * 30,
    age = 18 * 365,
    boosters = c(3, 12) * 30,
    booster_coverage = c(.9, .8),
    seasonal_boosters = TRUE
  )
  events <- create_events(parameters)

  booster_event <- mock_event(events$rtss_epi_boosters[[1]])

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

test_that('RTS,S EPI schedules for the following year with seasonal boosters', {
  timestep <- 5 * 365 + 4 * 30
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_rtss_epi(
    parameters,
    start = 10,
    end = timestep,
    coverage = 0.8,
    min_wait = 6 * 30,
    age = 18 * 365,
    boosters = c(3, 12) * 30,
    booster_coverage = c(.9, .8),
    seasonal_boosters = TRUE
  )
  events <- create_events(parameters)

  booster_event <- mock_event(events$rtss_epi_boosters[[1]])

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
