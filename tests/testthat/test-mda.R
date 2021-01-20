
test_that('MDA moves the diseased and non-diseased population correctly', {
  timestep <- 50
  renderer <- individual::Render$new(timestep)
  parameters <- get_parameters()
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

  variables <- list(
    state = mock_category(
      c('D', 'S', 'A', 'U', 'Tr'),
      c('D', 'S', 'A', 'U')
    ),
    birth = mock_double(-365 * c(2, 20, 5, 7)),
    infectivity = mock_double(c(.1, .2, .3, .4)),
    drug_time = mock_double(c(1, 2, 3, 4)),
    drug = mock_double(c(1, 2, 1, 2))
  )

  events$mda_administer <- mock_event(events$mda_administer)

  listener <- create_mda_listeners(
    variables,
    events$mda_administer,
    parameters$mda_drug,
    parameters$mda_end,
    parameters$mda_frequency,
    parameters$mda_min_age,
    parameters$mda_max_age,
    parameters$mda_coverage,
    get_correlation_parameters(parameters),
    'mda',
    parameters,
    renderer
  )

  mockery::stub(listener, 'bernoulli', mockery::mock(c(TRUE, TRUE)))
  mock_correlation <- mockery::mock(c(TRUE, TRUE))
  mockery::stub(listener, 'sample_intervention', mock_correlation)
  listener(timestep)

  expect_equal(
    mockery::mock_args(mock_correlation)[[1]][[1]],
    c(3, 4)
  )

  expect_bitset_update(
    variables$state$queue_update,
    'Tr',
    3
  )

  expect_bitset_update(
    variables$state$queue_update,
    'S',
    4,
    call = 2
  )

  expect_bitset_update(
    variables$infectivity$queue_update,
    c(.3, .4) * SP_AQ_params[[2]],
    c(3, 4)
  )

  expect_bitset_update(
    variables$drug$queue_update,
    1,
    c(3, 4)
  )

  expect_bitset_update(
    variables$drug_time$queue_update,
    50,
    c(3, 4)
  )

  mockery::expect_args(
    events$mda_administer$schedule,
    1,
    100
  )
})
