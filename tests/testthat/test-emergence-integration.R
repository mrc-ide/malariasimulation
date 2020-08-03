test_that('emergence process fails when there are not enough individuals', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  emergence_process <- create_mosquito_emergence_process(
    individuals$mosquito,
    list(list(), list()),
    states$Unborn,
    states$Sm,
    variables$mosquito_variety,
    2
  )

  api <- mock_api(
    list(
      mosquito = list(
        Sm = seq_len(1000),
        Im = seq_len(1000) + 1000,
        Unborn = c()
      )
    )
  )
  expect_error(
    with_mock(
      'malariasimulation:::mosquito_model_get_states' = mockery::mock(
        c(1000, 500, 100),
        c(1000, 500, 100)
      ),
      emergence_process(api)
    ),
    '*'
  )

  api <- mock_api(
    list(
      mosquito = list(
        Sm = seq_len(1000),
        Im = seq_len(1000) + 1000,
        Unborn = seq_len(1000) + 2000
      )
    )
  )
  expect_error(
    with_mock(
      'malariasimulation:::mosquito_model_get_states' = mockery::mock(
        c(100000, 50000, 10000),
        c(100000, 50000, 10000)
      ),
      emergence_process(api)
    ),
    '*'
  )
})

test_that('emergence_process creates the correct number of susceptables', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  emergence_process <- create_mosquito_emergence_process(
    individuals$mosquito,
    list(list(), list()),
    states$Unborn,
    states$Sm,
    variables$variety,
    2
  )

  api <- mock_api(
    list(
      mosquito = list(
        Sm = seq_len(1000),
        Im = seq_len(1000) + 1000,
        Unborn = seq_len(100000) + 2000
      )
    ),
    parameters = parameters
  )

  with_mock(
    'malariasimulation:::mosquito_model_get_states' = mockery::mock(
      c(100000, 50000, 10000),
      c(100000, 50000, 10000)
    ),
    emergence_process(api)
  )
  emergence_rate <- .5 * (1. - exp(-1./2))
  expected_target <- seq(emergence_rate * 20000) + 2000
  expected_variety <- c(
    rep(1, emergence_rate * 10000),
    rep(2, emergence_rate * 10000)
  )
  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$mosquito,
    states$Sm,
    expected_target
  )
  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$mosquito,
    variables$variety,
    expected_variety,
    expected_target
  )
})
