test_that('mortality_process resets humans correctly', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  mortality_process <- create_mortality_process(
    individuals$human,
    states$D,
    variables,
    events
  )

  api <- mock_api(
    list(
      human = list(
        D = c(1, 2),
        age = c(20, 24, 5, 39),
        is_severe = c(1., 1., 0., 0.),
        xi_group = c(1, 1, 2, 2),
        ICM = c(1, 2, 3, 4),
        IVM = c(1, 2, 3, 4)
      )
    ),
    timestep = 2
  )

  # NOTE: `with_mock` preferred here as `stub` suffers from locked binding issues
  with_mock(
    sample = mockery::mock(c(1), c(4)),
    'malariasimulation:::bernoulli' = mockery::mock(
      c(FALSE, FALSE, FALSE, TRUE),
      c(FALSE, TRUE)
    ),
    mortality_process(api)
  )

  died <- c(2, 4)

  expect(is.numeric(parameters$pvm), 'Pvm is not set')
  expect(is.numeric(parameters$pcm), 'Pcm is not set')

  update_args <- mockery::mock_args(api$queue_variable_update)
  cleared_args <- mockery::mock_args(api$clear_schedule)

  for (update in update_args) {
    expect_equal(update[[1]]$name, 'human')
    expect_setequal(update[[4]], died)
  }

  expect_setequal(
    vapply(update_args, function(update) update[[2]]$name, character(1)),
    c(
      'birth',
      'last_bitten',
      'last_infected',
      'ICM',
      'IVM',
      'IB',
      'ICA',
      'IVA',
      'ID',
      'is_severe'
    )
  )

  expect_setequal(
    vapply(cleared_args, function(cleared) cleared[[1]]$name, character(1)),
    c(
      'infection',
      'asymptomatic_infection',
      'subpatent_infection',
      'recovery'
    )
  )

})
