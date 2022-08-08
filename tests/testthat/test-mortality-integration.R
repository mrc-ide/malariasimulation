test_that('mortality_process resets humans correctly', {
  timestep <- 2
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_drugs(parameters, list(SP_AQ_params))
  parameters <- set_mda(
    parameters,
    drug = 1,
    timesteps = c(50, 150),
    coverages= c(1, 1),
    min_ages = c(5 * 365, 5 * 365),
    max_ages = c(10 * 365, 10 * 365)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  variables$state <- mock_category(
    c('S', 'A', 'D', 'U', 'Tr'),
    c('D', 'Tr', 'S', 'S')
  )
  variables$zeta_group = individual::CategoricalVariable$new(
    c('1', '2', '3', '4', '5'),
    c('1', '1', '2', '2')
  )
  variables$birth <- mock_double(-c(20, 24, 5, 39) * 365)
  variables$icm <- mock_double(c(1, 2, 3, 4))
  variables$ivm <- mock_double(c(1, 2, 3, 4))
  variables$ica <- mock_double(c(1, 2, 3, 4))
  variables$iva <- mock_double(c(1, 2, 3, 4))
  renderer <- individual::Render$new(timestep)

  mortality_process <- create_mortality_process(
    variables,
    events,
    renderer,
    parameters
  )

  mockery::stub( # natural deaths
    mortality_process,
    'bernoulli',
    c(2, 4)
  )

  mortality_process(timestep)

  expect_bitset_update(variables$state$queue_update_mock(), 'S', c(2, 4))
})

test_that('mortality_process samples deaths from a custom demography', {
  timestep <- 2
  parameters <- get_parameters()
  parameters$human_population <- 4
  ages <- c(50, 100) * 365
  deaths <- c(.5, .75)
  parameters <- set_demography(
    parameters,
    agegroups = ages,
    timesteps = 0,
    deathrates = matrix(deaths, nrow=1, ncol=2)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  variables$state <- mock_category(
    c('S', 'A', 'D', 'U', 'Tr'),
    c('D', 'Tr', 'S', 'S')
  )
  variables$zeta_group = individual::CategoricalVariable$new(
    c('1', '2', '3', '4', '5'),
    c('1', '1', '2', '2')
  )
  variables$birth <- mock_double(-c(20, 24, 75, 101) * 365)
  variables$icm <- mock_double(c(1, 2, 3, 4))
  variables$ivm <- mock_double(c(1, 2, 3, 4))
  variables$ica <- mock_double(c(1, 2, 3, 4))
  variables$iva <- mock_double(c(1, 2, 3, 4))
  renderer <- individual::Render$new(timestep)

  mortality_process <- create_mortality_process(
    variables,
    events,
    renderer,
    parameters
  )

  mortality_rng <- mockery::mock(c(2, 4))

  mockery::stub( # natural deaths
    mortality_process,
    'bernoulli_multi_p',
    mortality_rng
  )

  mortality_process(timestep)

  mockery::expect_args(mortality_rng, 1, c(.5, .5, .75, 1))
  expect_bitset_update(variables$state$queue_update_mock(), 'S', c(2, 4))
})

test_that('maternal immunity is sampled correctly', {
  timestep <- 2
  parameters <- get_parameters(list(human_population = 4))
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  variables$state <- mock_category(
    c('S', 'A', 'D', 'U', 'Tr'),
    c('D', 'Tr', 'S', 'S')
  )
  variables$zeta_group = individual::CategoricalVariable$new(
    c('1', '2', '3', '4', '5'),
    c('1', '1', '2', '2')
  )
  variables$birth <- mock_double(-c(20, 24, 5, 39) * 365)
  variables$icm <- mock_double(c(1, 2, 3, 4))
  variables$ivm <- mock_double(c(1, 2, 3, 4))
  variables$ica <- mock_double(c(1, 2, 3, 4))
  variables$iva <- mock_double(c(1, 2, 3, 4))
  renderer <- individual::Render$new(timestep)

  mockery::stub( # mothers
    sample_maternal_immunity,
    'sample.int',
    mockery::mock(
      1, # for 2
      1 # for 4
    )
  )

  sample_maternal_immunity(
    variables,
    individual::Bitset$new(4)$insert(c(2, 4)),
    timestep,
    parameters
  )

  expect_bitset_update(variables$icm$queue_update_mock(), parameters$pcm, 2)
  expect_bitset_update(variables$ivm$queue_update_mock(), parameters$pvm, 2)
  expect_bitset_update(variables$icm$queue_update_mock(), parameters$pcm, 4, call = 2)
  expect_bitset_update(variables$ivm$queue_update_mock(), parameters$pvm, 4, call = 2)
})
