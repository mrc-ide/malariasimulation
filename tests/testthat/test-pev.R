test_that('Mass vaccination strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = 10,
    coverages = 0.8,
    min_wait = 0,
    min_ages = 5 * 30,
    max_ages = 17 * 30,
    booster_spacing = c(18, 36) * 30,
    booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  expect_equal(parameters$pev, TRUE)
  expect_equal(parameters$mass_pev_timesteps, 10)
  expect_equal(parameters$mass_pev_coverages, .8)
  expect_equal(parameters$mass_pev_min_ages, 5 * 30)
  expect_equal(parameters$mass_pev_max_ages, 17 * 30)
  expect_equal(parameters$mass_pev_booster_spacing, c(18, 36) * 30)
  expect_equal(parameters$pev_profiles, list(rtss_profile, rtss_booster_profile, rtss_booster_profile))
  expect_equal(parameters$mass_pev_profile_indices, seq(3))

  expect_error(
    parameters <- set_mass_pev(
      parameters,
      profile = rtss_profile,
      timesteps = 10,
      coverages = -1, # less than 0
      min_wait = 0,
      min_ages = 5 * 30,
      max_ages = 17 * 30,
      booster_spacing = c(18, 36) * 30,
      booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    ), "all(coverages >= 0) && all(coverages <= 1) is not TRUE",
    fixed = TRUE
  )

  expect_error(
    parameters <- set_mass_pev(
      parameters,
      profile = rtss_profile,
      timesteps = 10,
      coverages = 1.5, # greater than 1
      min_wait = 0,
      min_ages = 5 * 30,
      max_ages = 17 * 30,
      booster_spacing = c(18, 36) * 30,
      booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    ), "all(coverages >= 0) && all(coverages <= 1) is not TRUE",
    fixed = TRUE
  )
})

test_that('set_mass_pev checks booster coverage matrix shape', {
  parameters <- get_parameters()
  expect_error(
    parameters <- set_mass_pev(
      parameters,
      profile = rtss_profile,
      coverages = c(0.1),
      timesteps = c(10),
      min_wait = 0,
      min_ages = 5 * 30,
      max_ages = 17 * 30,
      booster_spacing = c(18, 36) * 30,
      booster_coverage = matrix(c(.9, .8), nrow=2, ncol=1),
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    ),
    'booster_spacing, booster_coverage and booster_profile do not align',
    fixed = TRUE
  )
})

test_that('set_mass_pev checks booster_spacing is increasing', {
  parameters <- get_parameters()
  expect_error(
    parameters <- set_mass_pev(
      parameters,
      profile = rtss_profile,
      coverages = c(0.1),
      timesteps = c(10),
      min_wait = 0,
      min_ages = 5 * 30,
      max_ages = 17 * 30,
      booster_spacing = c(5, 5) * 30,
      booster_coverage = matrix(c(.9, .8), nrow=2, ncol=1),
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    ),
    'booster_spacing must be monotonically increasing',
    fixed = TRUE
  )
})

test_that('Mass vaccination fails pre-emptively for unaligned booster parameters', {
  parameters <- get_parameters()
  expect_error(
    set_mass_rtss(
      parameters,
      profile = rtss_profile,
      timesteps = 10,
      coverages = 0.8,
      min_wait = 0,
      min_ages = 5 * 30,
      max_ages = 17 * 30,
      booster_spacing = c(18, 36) * 30,
      booster_coverage = matrix(.9, nrow=1, ncol=1),
      booster_profile = list(rtss_booster_profile, rtss_booster_profile)
    )
  )
})

test_that('Infection considers pev efficacy', {
  timestep <- 100
  parameters <- get_parameters()
  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = 100,
    coverages = 0.8,
    min_wait = 0,
    min_ages = 5 * 30,
    max_ages = 17 * 30,
    booster_spacing = c(18, 36) * 30,
    booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  variables$state <- individual::CategoricalVariable$new(
    c('S', 'A', 'D', 'U', 'Tr'),
    rep('S', 4)
  )
  variables$birth <- individual::IntegerVariable$new(
    -c(8, 2.9, 3.2, 18.4) * 365 - 100
  )
  variables$last_eff_pev_timestep <- individual::IntegerVariable$new(
    c(-1, -1, 50, 50 + 30)
  )
  variables$pev_profile <- individual::IntegerVariable$new(
    c(-1, -1, 1, 2)
  )
  variables$drug <- individual::IntegerVariable$new(
    c(-1, -1, -1, -1)
  )
  variables$ib <- individual::IntegerVariable$new(
    rep(.2, 4)
  )

  infection_outcome <- CompetingOutcome$new(
    targeted_process = function(timestep, target){
      infection_process_resolved_hazard(timestep, target, variables, renderer, parameters)
    },
    size = parameters$human_population
  )

  # remove randomness from infection sampling
  bernoulli_mock <- mockery::mock(c(1, 2))
  mockery::stub(calculate_infections, 'bernoulli_multi_p', bernoulli_mock)

  # remove randomness from pev parameters
  mockery::stub(
    calculate_infections,
    'sample_pev_param',
    function(index, profiles, name) {
      vnapply(index, function(i) profiles[[i]][[name]][[1]]) # return mu
    },
    depth = 4
  )

  infection_rates <- calculate_infections(
    variables = variables,
    bitten_humans = individual::Bitset$new(4)$insert(seq(4)),
    parameters = parameters,
    renderer = mock_render(timestep),
    timestep = timestep,
    infection_outcome
  )

  expect_equal(
    rate_to_prob(infection_rates[infection_rates!=0]),
    c(0.590, 0.590, 0.215, 0.244),
    tolerance=1e-3
  )
})

test_that('Mass vaccinations update vaccination time', {
  timestep <- 100
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = c(100, 100 + 365),
    coverages = rep(0.8, 2),
    min_wait = 0,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    booster_spacing = c(18, 36) * 30,
    booster_coverage = matrix(c(.9, .8, .9, .8), nrow=2, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$birth <- individual::IntegerVariable$new(
    -c(18.3, 8, 2.9, 3.2, 18.4) * 365 + 100
  )
  variables$pev_timestep <- mock_integer(
    c(-1, -1, -1, 50, 50)
  )

  events$mass_pev_doses <- lapply(events$mass_pev_doses, mock_event)

  listener <- create_mass_pev_listener(
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
    events$mass_pev_doses[[1]]$schedule,
    1,
    c(1, 3),
    parameters$pev_doses[[1]]
  )

  mockery::expect_args(
    events$mass_pev_doses[[2]]$schedule,
    1,
    c(1, 3),
    parameters$pev_doses[[2]]
  )

  mockery::expect_args(
    events$mass_pev_doses[[3]]$schedule,
    1,
    c(1, 3),
    parameters$pev_doses[[3]]
  )
})

test_that('Mass vaccinations ignore EPI individuals', {
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
    min_wait = 2*365,
    age = 18 * 365,
    booster_spacing = c(18, 36) * 30,
    booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$birth <- individual::IntegerVariable$new(
    -c(18.3, 8, 2.9, 3.2, 18.4) * 365 + 100
  )
  variables$pev_timestep <- mock_integer(
    c(-1, -1, -1, 50, 50)
  )

  correlations <- get_correlation_parameters(parameters)

  listener <- create_mass_pev_listener(
    variables,
    events,
    parameters,
    correlations
  )

  sample_mock <- mockery::mock(c(TRUE, TRUE, FALSE, FALSE))

  mockery::stub(
    listener,
    'sample_intervention',
    sample_mock
  )
  
  # schedule id #1 for epi vaccination
  events$pev_epi_doses[[1]]$schedule(1, 0)

  listener(timestep)

  mockery::expect_args(
    sample_mock,
    1,
    c(3, 4, 5),
    'pev',
    .8,
    correlations
  )
})


test_that('Mass boosters update profile params and reschedule correctly', {
  parameters <- get_parameters()
  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = c(50, 50 + 365),
    coverages = rep(0.8, 2),
    min_wait = 0,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    booster_spacing = c(1, 6) * 30,
    booster_coverage = matrix(1, nrow=2, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  timestep <- 50 + 30

  variables$birth <- individual::IntegerVariable$new(
    -c(2.9, 3.2, 18.4) * 365 + 100
  )
  variables$last_pev_timestep <- mock_double(
    c(50, 50, 50)
  )
  variables$last_eff_pev_timestep <- mock_double(
    c(50, 50, 50)
  )
  variables$pev_profile <- mock_integer(
    c(1, 1, 1)
  )
  events$mass_pev_boosters <- lapply(events$mass_pev_boosters, mock_event)

  listener <- create_pev_booster_listener(
    variables = variables,
    coverage = parameters$mass_pev_booster_coverage,
    parameters$mass_pev_timesteps,
    booster_number = 1,
    pev_profile_index = 2,
    next_booster_event = events$mass_pev_boosters[[2]],
    next_booster_delay = 5 * 30,
    renderer = mock_render(timestep),
    strategy = 'mass'
  )

  mockery::stub(
    listener,
    'sample_bitset',
    mockery::mock(individual::Bitset$new(3)$insert(c(1, 2, 3)))
  )
  listener(timestep, individual::Bitset$new(3)$insert(c(1, 2, 3)))

  expect_bitset_update(
    variables$last_pev_timestep$queue_update_mock(),
    timestep,
    c(1, 2, 3)
  )

  expect_bitset_update(
    variables$last_eff_pev_timestep$queue_update_mock(),
    timestep,
    c(1, 2, 3)
  )

  expect_bitset_update(
    variables$pev_profile$queue_update_mock(),
    2,
    c(1, 2, 3)
  )

  expect_bitset_schedule(
    events$mass_pev_boosters[[2]]$schedule,
    c(1, 2, 3),
    5 * 30
  )
})

test_that('Mass booster coverages sample subpopulations correctly', {
  parameters <- get_parameters()
  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = 50,
    coverages = 0.8,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    min_wait = 0,
    booster_spacing = c(1, 6) * 30,
    booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )

  timestep <- 50 + 30

  events <- create_events(parameters)
  events$mass_pev_boosters <- lapply(events$mass_pev_boosters, mock_event)
  variables <- create_variables(parameters)

  variables$birth <- individual::IntegerVariable$new(
    -c(2.9, 3.2, 18.4) * 365 + 100
  )
  variables$last_pev_timestep <- mock_double(
    c(50, 50, 50)
  )
  variables$last_eff_pev_timestep <- mock_double(
    c(50, 50, 50)
  )
  variables$pev_profile <- mock_integer(
    c(1, 1, 1)
  )

  listener <- create_pev_booster_listener(
    variables = variables,
    coverage = parameters$mass_pev_booster_coverage,
    pev_distribution_timesteps = parameters$mass_pev_timesteps,
    booster_number = 1,
    pev_profile_index = 2,
    next_booster_event = events$mass_pev_boosters[[2]],
    next_booster_delay = 5 * 30,
    renderer = mock_render(timestep),
    strategy = 'mass'
  )

  sample_mock <- mockery::mock(individual::Bitset$new(3)$insert(c(2, 3)))

  mockery::stub(listener, 'sample_bitset', sample_mock)
  target <- individual::Bitset$new(3)$insert(c(1, 2, 3))
  listener(timestep, target)

  mockery::expect_args(sample_mock, 1, target, .9)

  expect_bitset_update(
    variables$last_pev_timestep$queue_update_mock(),
    timestep,
    c(2, 3)
  )

  expect_bitset_update(
    variables$last_eff_pev_timestep$queue_update_mock(),
    timestep,
    c(2, 3)
  )

  expect_bitset_update(
    variables$pev_profile$queue_update_mock(),
    2,
    c(2, 3)
  )

  expect_bitset_schedule(
    events$mass_pev_boosters[[2]]$schedule,
    c(2, 3),
    5 * 30
  )
})

test_that('mass pev targets correct age and respects min_wait', {
  timestep <- 5*365 
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = c(4, 5) * 365,
    coverages = c(0.8, 0.8),
    min_ages = 0,
    max_ages = 19 * 365,
    min_wait = 2*365,
    booster_spacing = c(1, 6) * 30,
    booster_coverage = matrix(c(.9, .8, .9, .8), nrow=2, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$birth <- individual::IntegerVariable$new(
    -c(18, 18, 30, 18, 18) * 365 + timestep
  )
  variables$last_pev_timestep <- mock_integer(
    c(50, -1, -1, 4*365, -1)
  )

  variables$pev_profile <- mock_integer(
    c(1, -1, -1, 1, -1)
  )

  correlations <- get_correlation_parameters(parameters)
  listener <- create_mass_pev_listener(
    variables,
    events,
    parameters,
    get_correlation_parameters(parameters)
  )

  sample_mock <- mockery::mock(c(TRUE, TRUE, FALSE))
  mockery::stub(
    listener,
    'sample_intervention',
    sample_mock
  )

  listener(timestep)

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

test_that('Mass efficacy listener works correctly', {
  timestep <- 50
  parameters <- get_parameters()
  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = 50,
    coverages = 0.8,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    min_wait = 0,
    booster_spacing = c(1, 6) * 30,
    booster_coverage = matrix(c(.9, .8), nrow=1, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile)
  )

  variables <- create_variables(parameters)
  variables$last_eff_pev_timestep <- mock_integer(c(-1, -1, -1))
  variables$pev_profile <- mock_integer(c(-1, -1, -1))
  listener <- create_pev_efficacy_listener(variables, 1)

  listener(timestep, individual::Bitset$new(3)$insert(c(1, 2, 3)))

  # vaccinated time
  expect_bitset_update(
    variables$last_eff_pev_timestep$queue_update_mock(),
    timestep,
    c(1, 2, 3)
  )

  # vaccinated profile
  expect_bitset_update(
    variables$pev_profile$queue_update_mock(),
    1,
    c(1, 2, 3)
  )
})

test_that('Mass dose events are not ruined by lazy evaluation', {
  parameters <- get_parameters(list(pev_doses = c(0, 7)))

  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = 50,
    coverages = 0.8,
    min_ages = c(1, 2, 3, 18) * 365,
    max_ages = (c(1, 2, 3, 18) + 1) * 365 - 1,
    min_wait = 0,
    booster_spacing = c(1, 6, 12) * 30,
    booster_coverage = matrix(c(.9, .8, .7), nrow=1, ncol=3),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile, rtss_booster_profile)
  )

  events <- create_events(parameters)
  variables <- create_variables(parameters)

  attach_pev_dose_listeners(
    variables = variables,
    parameters = parameters,
    pev_distribution_timesteps = parameters$mass_pev_timesteps,
    dose_events = events$mass_pev_doses,
    booster_events = events$mass_pev_boosters,
    booster_delays = parameters$mass_pev_booster_spacing,
    booster_coverages = parameters$mass_pev_booster_coverage,
    pev_profile_indices = parameters$mass_pev_profile_indices,
    strategy = 'mass',
    renderer = mock_render(1)
  )

  expect_equal(
    as.list(environment(
      events$mass_pev_boosters[[1]]$.listeners[[1]]
    ))$next_booster_event,
    events$mass_pev_boosters[[2]]
  )
  expect_equal(
    as.list(environment(
      events$mass_pev_boosters[[1]]$.listeners[[1]]
    ))$next_booster_delay,
    5 * 30
  )
  expect_equal(
    as.list(environment(
      events$mass_pev_boosters[[1]]$.listeners[[1]]
    ))$booster_number,
    1
  )
})

test_that('pev Antibodies are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_pev_antibodies(
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
    calculate_pev_efficacy(
      c(400, 400, 100, 50),
      rtss_profile$vmax,
      rtss_profile$beta,
      rtss_profile$alpha
    ),
    c(0.685, 0.685, 0.466, 0.349),
    tolerance = 1e-3
  )
})

test_that('pev timed booster coverage can select the first coverage for the first booster', {
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
  )
  events <- create_events(parameters)

  booster_event <- mock_event(events$pev_epi_boosters[[1]])

  listener <- create_pev_booster_listener(
    variables = create_variables(parameters),
    coverage = parameters$pev_epi_booster_coverage,
    pev_distribution_timesteps = parameters$pev_epi_timesteps,
    booster_number = 1,
    pev_profile_index = 2,
    next_booster_event = booster_event,
    next_booster_delay = 9 * 30,
    renderer = mock_render(timestep),
    strategy = 'epi'
  )

  target <- individual::Bitset$new(5)$insert(seq(5))

  mock_sample_bitset = mockery::mock(individual::Bitset$new(5)$insert(c(1, 2)))
  mockery::stub(
    listener,
    'sample_bitset',
    mock_sample_bitset
  )

  listener(timestep, target)

  mockery::expect_args(
    mock_sample_bitset,
    1,
    target,
    .9
  )
})


test_that('pev boosters can select the second coverage for the first booster', {
  timestep <- 5 * 365
  parameters <- get_parameters(list(human_population = 5))
  parameters <- set_pev_epi(
    parameters,
    profile = rtss_profile,
    timesteps = c(10, 30),
    coverages = c(0.8, 0.4),
    min_wait = 6 * 30,
    age = 18 * 365,
    booster_spacing = c(3, 12) * 30,
    booster_coverage = matrix(c(.9, .45, .8, .8), nrow=2, ncol=2),
    booster_profile = list(rtss_booster_profile, rtss_booster_profile),
  )
  events <- create_events(parameters)

  booster_event <- mock_event(events$pev_epi_boosters[[1]])

  listener <- create_pev_booster_listener(
    variables = create_variables(parameters),
    coverage = parameters$pev_epi_booster_coverage,
    pev_distribution_timesteps = parameters$pev_epi_timesteps,
    booster_number = 1,
    pev_profile_index = 2,
    next_booster_event = booster_event,
    next_booster_delay = 9 * 30,
    renderer = mock_render(timestep),
    strategy = 'epi'
  )

  target <- individual::Bitset$new(5)$insert(seq(5))

  mock_sample_bitset = mockery::mock(individual::Bitset$new(5)$insert(c(1, 2)))
  mockery::stub(
    listener,
    'sample_bitset',
    mock_sample_bitset
  )

  listener(timestep, target)

  mockery::expect_args(
    mock_sample_bitset,
    1,
    target,
    .45
  )
})
