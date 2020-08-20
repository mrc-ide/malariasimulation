test_that('mortality_process resets humans correctly', {
  parameters <- get_parameters(list(severe_enabled = 1))
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
  individuals <- create_individuals(
    states,
    variables,
    events,
    parameters
  )

  mortality_process <- create_mortality_process(
    individuals$human,
    states$D,
    states$Tr,
    variables,
    events
  )

  api <- mock_api(
    list(
      human = list(
        D = c(1),
        Tr = c(2),
        is_severe = c(1., 1., 0., 0.),
        zeta_group = c(1, 1, 2, 2),
        ICM = c(1, 2, 3, 4),
        IVM = c(1, 2, 3, 4)
      )
    ),
    parameters = parameters,
    timestep = 2
  )

  with_mock(
    sample = mockery::mock(c(1), c(4)),
    'malariasimulation:::bernoulli' = mockery::mock(
      c(4),
      numeric(0),
      c(1)
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
      'last_boosted_ib',
      'last_boosted_ica',
      'last_boosted_iva',
      'last_boosted_id',
      'ICM',
      'IVM',
      'IB',
      'ICA',
      'IVA',
      'ID',
      'drug',
      'drug_time',
      'infectivity'
    )
  )

  expect_setequal(
    vapply(cleared_args, function(cleared) cleared[[1]]$name, character(1)),
    c(
      'infection',
      'asymptomatic_infection',
      'mda_administer'
    )
  )
})
