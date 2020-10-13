
test_that('MDA moves the diseased and non-diseased population correctly', {
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
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  create_event_based_processes(
    individuals,
    states,
    variables,
    events,
    parameters
  )

  api <- mock_api(
    list(
      human = list(
        U = c(1),
        A = c(2),
        D = c(3),
        S = c(4),
        birth = -365 * c(2, 20, 5, 7),
        infectivity = c(.1, .2, .3, .4)
      )
    ),
    timestep = 50,
    parameters = parameters
  )

  listener <- events$mda_administer$listeners[[1]]
  mockery::stub(listener, 'bernoulli', mockery::mock(c(TRUE, TRUE)))
  mock_correlation <- mockery::mock(c(TRUE, TRUE))
  mockery::stub(listener, 'sample_intervention', mock_correlation)
  listener(api, c(1))

  expect_equal(
    mockery::mock_args(mock_correlation)[[1]][[1]],
    c(3, 4)
  )

  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$human,
    states$Tr,
    3
  )

  mockery::expect_args(
    api$queue_state_update,
    2,
    individuals$human,
    states$S,
    4
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$infectivity,
    c(.3, .4) * SP_AQ_params[[2]],
    c(3, 4)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$drug,
    1,
    c(3, 4)
  )

  mockery::expect_args(
    api$queue_variable_update,
    3,
    individuals$human,
    variables$drug_time,
    50,
    c(3, 4)
  )

  mockery::expect_args(
    api$schedule,
    1,
    events$mda_administer,
    c(1), #Dummy target
    100
  )
})

test_that('SMC moves the diseased and non-diseased population correctly', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(SP_AQ_params))
  parameters <- set_smc(
    parameters,
    1, # drug
    50, # start timestep
    50 + 365 * 5, # last timestep
    100, # frequency
    0, # min age
    6 * 365, # max age
    1 # coverage
  )
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  create_event_based_processes(
    individuals,
    states,
    variables,
    events,
    parameters
  )

  api <- mock_api(
    list(
      human = list(
        U = c(1),
        A = c(2),
        D = c(3),
        S = c(4),
        birth = c(2, -365 * 20, -365 * 5, -365 * 7),
        infectivity = c(1, 2, 3, 4)
      )
    ),
    timestep = 50,
    parameters = parameters
  )

  listener <- events$smc_administer$listeners[[1]]
  mockery::stub(listener, 'bernoulli', mockery::mock(c(TRUE, FALSE)))
  mock_correlation <- mockery::mock(c(TRUE, TRUE))
  mockery::stub(listener, 'sample_intervention', mock_correlation)
  listener(api, c(1))

  expect_equal(
    mockery::mock_args(mock_correlation)[[1]][[1]],
    c(1, 3)
  )

  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$human,
    states$S,
    1
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$infectivity,
    1 * SP_AQ_params[[2]],
    1 
  )

  mockery::expect_args(
    api$schedule,
    1,
    events$smc_administer,
    c(1),
    100
  )
})
