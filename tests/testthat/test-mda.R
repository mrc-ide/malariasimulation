
test_that('mda initialises with full coverage correctly', {
  parameters <- get_parameters()
  parameters <- add_drug(parameters, 1)
  parameters <- add_mda(
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
  mda_events <- create_mda_events(parameters)
  individuals <- create_individuals(states, variables, events, parameters, mda_events)
  listeners <- create_mda_listeners(
    individuals$human,
    states,
    variables,
    1,
    mda_events[[1]]
  )

  api <- mock_api(
    list(
      human = list(
        U = c(1),
        A = c(2),
        S = c(4),
        D = c(3),
        birth = c(2, -365 * 20, -365 * 5, -365 * 7)
      )
    ),
    timestep = 50,
    parameters = parameters
  )

  m <- mockery::mock()
  mockery::stub(listeners$mda_enrollment_listener, 'mda_administer_listener', m)
  listeners$mda_enrollment_listener(api, NULL)

  mockery::expect_args(
    m,
    1,
    api,
    c(3, 4),
    1
  )
})

test_that('mda moves the diseased and non-diseased population correctly', {
  parameters <- get_parameters()
  parameters <- add_drug(parameters, 1)
  parameters <- add_mda(
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
  mda_events <- create_mda_events(parameters)
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(
    states,
    variables,
    events,
    parameters,
    mda_events
  )
  listeners <- create_mda_listeners(
    individuals$human,
    states,
    variables,
    1,
    mda_events[[1]]
  )

  api <- mock_api(
    list(
      human = list(
        U = c(1),
        A = c(2),
        D = c(3),
        S = c(4),
        birth = c(2, -365 * 20, -365 * 5, -365 * 7)
      )
    ),
    timestep = 50,
    parameters = parameters
  )

  listeners$mda_administer_listener(api, c(3, 4))

  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$human,
    states$T,
    3
  )

  mockery::expect_args(
    api$queue_state_update,
    2,
    individuals$human,
    states$Ph,
    4
  )

  mockery::expect_args(
    api$schedule,
    1,
    mda_events[[1]]$mda_administer,
    c(3, 4),
    100
  )
})

test_that('multiple mdas can be scheduled', {
  parameters <- get_parameters()
  parameters <- add_drug(parameters, .7)
  parameters <- add_drug(parameters, .6)
  parameters <- add_mda(
    parameters,
    1, # drug
    50, # start timestep
    50 + 365 * 5, # last timestep
    100, # frequency
    5 * 365, # min age
    10 * 365, # max age
    1 # coverage
  )
  parameters <- add_mda(
    parameters,
    2, # drug
    50, # start timestep
    50 + 365 * 5, # last timestep
    200, # frequency
    10 * 365, # min age
    30 * 365, # max age
    1 # coverage
  )
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  mda_events <- create_mda_events(parameters)
  individuals <- create_individuals(
    states,
    variables,
    events,
    parameters,
    mda_events
  )
  listeners <- lapply(
    seq_along(mda_events),
    function(i) {
      create_mda_listeners(
        individuals$human,
        states,
        variables,
        i,
        mda_events[[i]]
      )
    }
  )

  api <- mock_api(
    list(
      human = list(
        U = c(1),
        A = c(2),
        D = c(3),
        S = c(4),
        birth = c(2, -365 * 20, -365 * 5, -365 * 7)
      )
    ),
    timestep = 50,
    parameters = parameters
  )

  bernoulli_mock <- mockery::mock(rep(TRUE, 2), rep(TRUE, 1))
  mockery::stub(
    listeners[[1]]$mda_administer_listener,
    'bernoulli',
    bernoulli_mock
  )
  mockery::stub(
    listeners[[2]]$mda_administer_listener,
    'bernoulli',
    bernoulli_mock
  )
  listeners[[1]]$mda_administer_listener(api, c(3, 4))
  listeners[[2]]$mda_administer_listener(api, c(2))

  mockery::expect_args(bernoulli_mock, 1, 2, .7)
  mockery::expect_args(bernoulli_mock, 2, 1, .6)

  mockery::expect_args(
    api$schedule,
    1,
    mda_events[[1]]$mda_administer,
    c(3, 4),
    100
  )

  mockery::expect_args(
    api$schedule,
    2,
    mda_events[[2]]$mda_administer,
    2,
    200
  )
})
