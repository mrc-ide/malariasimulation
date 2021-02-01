test_that('mortality_process resets humans correctly', {
  timestep <- 2
  parameters <- get_parameters(list(severe_enabled = 1, human_population = 4))
  parameters <- set_drugs(parameters, list(SP_AQ_params))
  parameters <- set_mda(
    parameters,
    1, # drug
    50, # start timestep
    50 + 365 * 5, # last timestep
    100, # frequency
    5 * 365, # min age
    10 * 365, # max age
    1 # coverage
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  variables$state <- individual::CategoricalVariable$new(
    c('S', 'A', 'D', 'U', 'Tr'),
    c('D', 'Tr', 'S', 'S')
  )
  variables$is_severe = individual::CategoricalVariable$new(
    c('yes', 'no'),
    c('yes', 'yes', 'no', 'no')
  )
  variables$zeta_group = individual::CategoricalVariable$new(
    c('1', '2', '3', '4', '5'),
    c('1', '1', '2', '2')
  )
  variables$birth <- mock_double(-c(20, 24, 5, 39) * 365)
  variables$icm <- mock_double(c(1, 2, 3, 4))
  variables$ivm <- mock_double(c(1, 2, 3, 4))
  renderer <- individual::Render$new(timestep)

  mortality_process <- create_mortality_process(
    variables,
    events,
    renderer,
    parameters
  )

  mockery::stub(mortality_process, 'sample', mockery::mock(1, 4), depth = 2)
  first_sample <- individual::Bitset$new(4)
  first_sample$insert(4)
  second_sample <- individual::Bitset$new(4)
  third_sample <- individual::Bitset$new(4)
  third_sample$insert(1)
  mockery::stub(
    mortality_process,
    'sample_bitset',
    mockery::mock(first_sample, second_sample, third_sample),
    depth = 2
  )

  mortality_process(timestep)

  died <- c(2, 4)

  expect_bitset_update(variables$state, 'S', died)
  expect_bitset_update(variables$is_severe, 'no', died)
})
