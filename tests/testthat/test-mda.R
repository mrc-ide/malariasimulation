test_that('MDA inputs align', {
  parameters <- get_parameters()
  expect_error(
    parameters <- set_mda(
      parameters,
      drug = 1,
      timesteps = c(50, 150),
      coverages= 1,
      min_ages = c(5* 365, 5 * 365),
      max_ages = c(10 * 365, 10 * 365)
    ), "coverages and timesteps do no align")

  expect_error(
    parameters <- set_mda(
      parameters,
      drug = 1,
      timesteps = c(50, 150),
      coverages = c(1, 1),
      min_ages = 5 * 365,
      max_ages = c(10 * 365, 10 * 365)
    ), "minimum ages and timesteps do no align")

  expect_error(
    parameters <- set_mda(
      parameters,
      drug = 1,
      timesteps = c(50, 150),
      coverages= c(1, 1),
      min_ages = c(5* 365, 5 * 365),
      max_ages = 10 * 365
    ), "maximum ages and timesteps do no align")

  expect_error(
    parameters <- set_mda(
      parameters,
      drug = 1,
      timesteps = c(50, 150),
      coverages= c(-1, 0.5),
      min_ages = c(5* 365, 5 * 365),
      max_ages = c(10 * 365, 10 * 365)
    ), "all(coverages >= 0) && all(coverages <= 1) is not TRUE",
    fixed = TRUE)

  expect_error(
    parameters <- set_mda(
      parameters,
      drug = 1,
      timesteps = c(50, 150),
      coverages= c(0.5, 1.5),
      min_ages = c(5* 365, 5 * 365),
      max_ages = c(10 * 365, 10 * 365)
    ), "all(coverages >= 0) && all(coverages <= 1) is not TRUE",
    fixed = TRUE)

})

test_that('MDA moves the diseased and non-diseased population correctly', {
  timestep <- 50
  renderer <- individual::Render$new(timestep)
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_drugs(parameters, list(SP_AQ_params_falciparum))
  parameters <- set_mda(
    parameters,
    drug = 1,
    timesteps = c(50, 150),
    coverages= c(1, 1),
    min_ages = c(5 * 365, 0),
    max_ages = c(10 * 365, 100 * 365)
  )
  events <- create_events(parameters)

  variables <- list(
    state = mock_category(
      c('D', 'S', 'A', 'U', 'Tr'),
      c('D', 'S', 'A', 'U')
    ),
    birth = mock_double(-365 * c(2, 20, 5, 7)),
    infectivity = mock_double(c(.1, .2, .3, .4)),
    id = mock_double(c(.1, .2, .3, .4)),
    drug_time = mock_double(c(1, 2, 3, 4)),
    drug = mock_double(c(1, 2, 1, 2))
  )

  events$mda_administer <- mock_event(events$mda_administer)

  listener <- create_mda_listeners(
    variables,
    parameters$mda_drug,
    parameters$mda_timesteps,
    parameters$mda_coverages,
    parameters$mda_min_ages,
    parameters$mda_max_ages,
    get_correlation_parameters(parameters),
    'mda',
    parameters,
    renderer
  )

  mockery::stub(listener, 'bernoulli', mockery::mock(c(TRUE, TRUE)))
  mock_correlation <- mockery::mock(c(TRUE, TRUE))
  mockery::stub(listener, 'sample_intervention', mock_correlation)
  listener(timestep)
  mockery::stub(
    listener,
    'calculate_asymptomatic_detectable',
    mockery::mock(individual::Bitset$new(4)$insert(3))
  )

  expect_equal(
    mockery::mock_args(mock_correlation)[[1]][[1]],
    c(3, 4)
  )
  expect_bitset_update(
    variables$state$queue_update_mock(),
    'Tr',
    3
  )

  expect_bitset_update(
    variables$state$queue_update_mock(),
    'S',
    4,
    call = 2
  )

  expect_bitset_update(
    variables$infectivity$queue_update_mock(),
    c(.3, .4) * SP_AQ_params_falciparum[2],
    c(3, 4)
  )

  expect_bitset_update(
    variables$drug$queue_update_mock(),
    1,
    c(3, 4)
  )

  expect_bitset_update(
    variables$drug_time$queue_update_mock(),
    50,
    c(3, 4)
  )

  mockery::expect_args(
    events$mda_administer$schedule,
    1,
    100
  )

})

test_that('MDA moves the diseased and non-diseased population correctly - second round, varying age range', {
  timestep <- 150
  renderer <- individual::Render$new(timestep)
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_drugs(parameters, list(SP_AQ_params_falciparum))
  parameters <- set_mda(
    parameters,
    drug = 1,
    timesteps = c(50, 150),
    coverages= c(1, 1),
    min_ages = c(5 * 365, 0),
    max_ages = c(10 * 365, 100 * 365)
  )
  events <- create_events(parameters)

  variables <- list(
    state = mock_category(
      c('D', 'S', 'A', 'U', 'Tr'),
      c('D', 'S', 'A', 'U')
    ),
    birth = mock_double(-365 * c(2, 20, 5, 7)),
    infectivity = mock_double(c(.1, .2, .3, .4)),
    id = mock_double(c(.1, .2, .3, .4)),
    drug_time = mock_double(c(1, 2, 3, 4)),
    drug = mock_double(c(1, 2, 1, 2))
  )

  events$mda_administer <- mock_event(events$mda_administer)

  listener <- create_mda_listeners(
    variables,
    parameters$mda_drug,
    parameters$mda_timesteps,
    parameters$mda_coverages,
    parameters$mda_min_ages,
    parameters$mda_max_ages,
    get_correlation_parameters(parameters),
    'mda',
    parameters,
    renderer
  )

  mockery::stub(listener, 'bernoulli', mockery::mock(c(TRUE, TRUE, TRUE, TRUE)))
  mock_correlation <- mockery::mock(c(TRUE, TRUE, TRUE, TRUE))
  mockery::stub(listener, 'sample_intervention', mock_correlation)
  mockery::stub(
    listener,
    'calculate_asymptomatic_detectable',
    mockery::mock(individual::Bitset$new(4)$insert(3))
  )
  listener(timestep)

  expect_equal(
    mockery::mock_args(mock_correlation)[[1]][[1]],
    c(1, 2, 3, 4)
  )
  expect_bitset_update(
    variables$state$queue_update_mock(),
    'Tr',
    c(1, 3)
  )

  expect_bitset_update(
    variables$state$queue_update_mock(),
    'S',
    c(2, 4),
    call = 2
  )

  expect_bitset_update(
    variables$infectivity$queue_update_mock(),
    c(.1, .2, .3, .4) * SP_AQ_params_falciparum[2],
    c(1, 2, 3, 4)
  )

  expect_bitset_update(
    variables$drug$queue_update_mock(),
    1,
    c(1, 2, 3, 4)
  )

  expect_bitset_update(
    variables$drug_time$queue_update_mock(),
    150,
    c(1, 2, 3, 4)
  )
})

test_that('MDA ignores non-detectable asymptomatics', {
  timestep <- 150
  renderer <- individual::Render$new(timestep)
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_drugs(parameters, list(SP_AQ_params_falciparum))
  parameters <- set_mda(
    parameters,
    drug = 1,
    timesteps = c(50, 150),
    coverages= c(1, 1),
    min_ages = c(5 * 365, 0),
    max_ages = c(10 * 365, 100 * 365)
  )
  events <- create_events(parameters)

  variables <- list(
    state = mock_category(
      c('D', 'S', 'A', 'U', 'Tr'),
      c('D', 'S', 'A', 'U')
    ),
    birth = mock_double(-365 * c(2, 20, 5, 7)),
    infectivity = mock_double(c(.1, .2, .3, .4)),
    id = mock_double(c(.1, .2, .3, .4)),
    drug_time = mock_double(c(1, 2, 3, 4)),
    drug = mock_double(c(1, 2, 1, 2))
  )

  events$mda_administer <- mock_event(events$mda_administer)

  listener <- create_mda_listeners(
    variables,
    parameters$mda_drug,
    parameters$mda_timesteps,
    parameters$mda_coverages,
    parameters$mda_min_ages,
    parameters$mda_max_ages,
    get_correlation_parameters(parameters),
    'mda',
    parameters,
    renderer
  )

  mockery::stub(listener, 'bernoulli', mockery::mock(c(TRUE, TRUE, TRUE, TRUE)))
  mock_correlation <- mockery::mock(c(TRUE, TRUE, TRUE, TRUE))
  mockery::stub(listener, 'sample_intervention', mock_correlation)
  mockery::stub(
    listener,
    'calculate_asymptomatic_detectable',
    mockery::mock(individual::Bitset$new(4))
  )
  listener(timestep)

  expect_bitset_update(
    variables$state$queue_update_mock(),
    'Tr',
    c(1)
  )

  expect_bitset_update(
    variables$state$queue_update_mock(),
    'S',
    c(2, 3, 4),
    call = 2
  )

  expect_bitset_update(
    variables$infectivity$queue_update_mock(),
    c(.1, .2, .3, .4) * SP_AQ_params_falciparum[2],
    c(1, 2, 3, 4)
  )

  expect_bitset_update(
    variables$drug$queue_update_mock(),
    1,
    c(1, 2, 3, 4)
  )

  expect_bitset_update(
    variables$drug_time$queue_update_mock(),
    150,
    c(1, 2, 3, 4)
  )
})

test_that('calculate_asymptomatic_detectable returns only detectable A individuals', {
  timestep <- 50
  state <- individual::CategoricalVariable$new(c('S', 'A'), c('S', 'A', 'A', 'A'))
  birth <- individual::IntegerVariable$new(timestep - c(0, 5, 30, 6) * 365)
  immunity <- individual::DoubleVariable$new(c(2.4, 1.2, 0., 4.))
  parameters <- get_parameters()

  bernoulli_mock <- mockery::mock(c(2, 3))

  mockery::stub(
    calculate_asymptomatic_detectable,
    'bernoulli_multi_p',
    bernoulli_mock
  )
  expect_equal(
    calculate_asymptomatic_detectable(
      state,
      birth,
      immunity,
      parameters,
      timestep
    )$to_vector(),
    c(3, 4)
  )
})
