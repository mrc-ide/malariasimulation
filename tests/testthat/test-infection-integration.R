
test_that('human infection_process works for non-severe clinical cases', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  infection_process <- create_infection_process(
    individuals,
    states,
    variables,
    events
  )

  parameters$severe_enabled = 1

  api <- mock_api(
    list(
      human = list(
        S = c(2),
        U = c(3),
        A = c(4),
        D = c(1),
        birth = -c(20, 24, 5, 39) * 365 + 5,
        IB = c(.2, .3, .5, .9),
        zeta = c(.2, .3, .5, .9),
        ICA = c(.2, .3, .5, .9),
        IVA = c(.2, .3, .5, .9),
        ICM = c(.2, .3, .5, .9),
        IVM = c(.2, .3, .5, .9),
        infectivity = c(.6, 0, .2, .3),
        drug = c(0, 0, 0, 0),
        last_boosted_ib = c(-1, -1, 1, -1),
        last_boosted_ica = c(-1, -1, 1, -1),
        last_boosted_iva = c(-1, -1, 1, -1),
        last_boosted_id = c(-1, -1, 1, -1),
        ID = c(.2, .3, .5, .9)
      ),
      mosquito = list(
        Im = 1:100,
        variety = c(rep(1, 25), rep(2, 25), rep(3, 50))
      )
    ),
    timestep = 5,
    parameters = parameters
  )

  bernoulli_mock <- mockery::mock(
    c(TRUE, FALSE, TRUE, FALSE), # bitten
    c(TRUE),                     # infected
    c(TRUE),                     # clinical
    c(FALSE)                     # severe
  )

  api$get_scheduled = mockery::mock(4, 1)

  with_mock(
    'malariasimulation:::bernoulli_multi_p' = bernoulli_mock,
    'malariasimulation:::bernoulli' = mockery::mock(numeric(0)), # mock seek treatment
    infection_process(api)
  )

  mockery::expect_called(bernoulli_mock, 4)

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$ib,
    1.2,
    1
  )
  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$last_boosted_ib,
    5,
    1
  )
  mockery::expect_args(
    api$queue_variable_update,
    3,
    individuals$human,
    variables$is_severe,
    FALSE,
    3
  )
  mockery::expect_args(
    api$schedule,
    1,
    events$infection,
    3,
    12
  )
})

test_that('human infection_process works for one drug', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  parameters <- set_clinical_treatment(parameters, .5, 2, 1)
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  infection_process <- create_infection_process(
    individuals,
    states,
    variables,
    events
  )

  parameters$severe_enabled = 1

  api <- mock_api(
    list(
      human = list(
        S = c(2),
        U = c(3),
        A = c(4),
        D = c(1),
        birth = -c(20, 24, 5, 39) * 365 + 5,
        IB = c(.2, .3, .5, .9),
        xi = c(.2, .3, .5, .9),
        ICA = c(.2, .3, .5, .9),
        IVA = c(.2, .3, .5, .9),
        ICM = c(.2, .3, .5, .9),
        IVM = c(.2, .3, .5, .9),
        infectivity = c(.6, 0, .2, .3),
        drug = c(0, 0, 0, 0),
        last_boosted_ib = c(-1, -1, 1, -1),
        last_boosted_ica = c(-1, -1, 1, -1),
        last_boosted_iva = c(-1, -1, 1, -1),
        last_boosted_id = c(-1, -1, 1, -1),
        ID = c(.2, .3, .5, .9)
      ),
      mosquito = list(
        Im = 1:100,
        variety = c(rep(1, 25), rep(2, 25), rep(3, 50))
      )
    ),
    timestep = 5,
    parameters = parameters
  )

  bernoulli_mock <- mockery::mock(
    c(FALSE, TRUE, TRUE, FALSE), # bitten
    c(TRUE, TRUE),               # infected
    c(TRUE, TRUE),               # clinical
    c(FALSE, TRUE),              # severe
    c(TRUE, TRUE)                # treatment successful
  )

  with_mock(
    'malariasimulation:::bernoulli_multi_p' = bernoulli_mock,
    'malariasimulation:::bernoulli' = mockery::mock(c(1, 2)),
    infection_process(api)
  )

  api$get_scheduled = mockery::mock(4, 1)

  infection_process(api)

  # expect treated individual to go to Tr state
  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$human,
    states$Tr,
    c(2, 3)
  )

  # expect rel_c to be applied correctly
  mockery::expect_args(
    api$queue_variable_update,
    10,
    individuals$human,
    variables$infectivity,
    rep(parameters$cd * 0.09434, 2),
    c(2, 3)
  )
  
  # expect drug and time to be applied correctly
  mockery::expect_args(
    api$queue_variable_update,
    11,
    individuals$human,
    variables$drug,
    c(2, 2),
    c(2, 3)
  )
  mockery::expect_args(
    api$queue_variable_update,
    12,
    individuals$human,
    variables$drug_time,
    5,
    c(2, 3)
  )
})

test_that('human infection_process works for two drugs', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  parameters <- set_clinical_treatment(parameters, .5, c(1, 2), c(.5, .5))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  infection_process <- create_infection_process(
    individuals,
    states,
    variables,
    events
  )

  parameters$severe_enabled = 1

  api <- mock_api(
    list(
      human = list(
        S = c(2),
        U = c(3),
        A = c(4),
        D = c(1),
        birth = -c(20, 24, 5, 39) * 365 + 5,
        IB = c(.2, .3, .5, .9),
        xi = c(.2, .3, .5, .9),
        ICA = c(.2, .3, .5, .9),
        IVA = c(.2, .3, .5, .9),
        ICM = c(.2, .3, .5, .9),
        IVM = c(.2, .3, .5, .9),
        infectivity = c(.6, 0, .2, .3),
        drug = c(0, 0, 0, 0),
        last_boosted_ib = c(-1, -1, 1, -1),
        last_boosted_ica = c(-1, -1, 1, -1),
        last_boosted_iva = c(-1, -1, 1, -1),
        last_boosted_id = c(-1, -1, 1, -1),
        ID = c(.2, .3, .5, .9)
      ),
      mosquito = list(
        Im = 1:100,
        variety = c(rep(1, 25), rep(2, 25), rep(3, 50))
      )
    ),
    timestep = 5,
    parameters = parameters
  )

  bernoulli_mock <- mockery::mock(
    c(TRUE, FALSE, TRUE, FALSE), # bitten
    c(TRUE),                     # infected
    c(TRUE),                     # clinical
    c(FALSE),                    # severe
    c(TRUE)                      # treatment successful
  )

  with_mock(
    'malariasimulation:::bernoulli_multi_p' = bernoulli_mock,
    'malariasimulation:::bernoulli' = mockery::mock(1),
    sample.int = mockery::mock(2), # Return DCH drug
    infection_process(api)
  )

  api$get_scheduled = mockery::mock(4, 1)

  infection_process(api)

  # expect treated individual to go to Tr state
  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$human,
    states$Tr,
    3
  )

  # expect rel_c to be applied correctly
  mockery::expect_args(
    api$queue_variable_update,
    4,
    individuals$human,
    variables$infectivity,
    parameters$cd * 0.09434,
    3
  )
  
  # expect drug and time to be applied correctly
  mockery::expect_args(
    api$queue_variable_update,
    5,
    individuals$human,
    variables$drug,
    2,
    3
  )
  mockery::expect_args(
    api$queue_variable_update,
    6,
    individuals$human,
    variables$drug_time,
    5,
    3
  )
})

test_that('prophylaxis is considered for medicated humans', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  ib = c(.2, .3, .5, .9)
  api <- mock_api(
    list(
      human = list(
        D = c(1),
        S = c(2),
        U = c(3),
        A = c(4),
        drug = c(0, 2, 1, 0),
        drug_time = c(-1, 49, 40, -1)
      )
    ),
    timestep = 50,
    parameters = parameters
  )
  m <- mockery::mock(c(TRUE, TRUE, TRUE, TRUE))

  with_mock(
    'malariasimulation:::bernoulli_multi_p' = m,
    calculate_infections(api, individuals$human, states, variables, seq(4), ib)
  )

  expect_equal(
    mockery::mock_args(m)[[1]][[2]],
    c(0.590, 0.590, 0.384),
    tolerance = 1e-3
  )
})

test_that('boost_immunity respects the delay period', {
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(
    states,
    variables,
    create_events(),
    parameters
  )

  level <- c(2.4, 1.2, 0., 4.)
  last_boosted <- c(11, 5, 1, 13)
  index <- seq(4)
  timestep <- 15
  delay <- 4

  api <- mock_api(
    list(
      human = list(
        ID = level,
        last_boosted_id = last_boosted
      )
    )
  )

  boost_immunity(
    api,
    individuals$human,
    variables$id,
    index,
    level,
    variables$last_boosted_id,
    timestep,
    delay
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$id,
    c(3.4, 2.2, 1),
    seq(3)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$last_boosted_id,
    15,
    seq(3)
  )
})

test_that('boost_immunity respects the delay period', {
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(
    states,
    variables,
    create_events(),
    parameters
  )

  level <- c(2.4, 1.2, 0., 4., 0.)
  last_boosted <- c(11, 5, 1, 13, -1)
  index <- seq(5)
  timestep <- 15
  delay <- 4

  api <- mock_api(
    list(
      human = list(
        ID = level,
        last_boosted_id = last_boosted
      )
    )
  )

  boost_immunity(
    api,
    individuals$human,
    variables$id,
    index,
    level,
    variables$last_boosted_id,
    timestep,
    delay
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$id,
    c(3.4, 2.2, 1, 1),
    c(seq(3), 5)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$last_boosted_id,
    15,
    c(seq(3), 5)
  )
})

test_that('boost_immunity does not update when there is no-one to update', {
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(
    states,
    variables,
    create_events(),
    parameters
  )

  level <- c(2.4, 1.2, 0., 4., 0.)
  last_boosted <- c(12, 14, 14, 13, 13)
  index <- seq(5)
  timestep <- 15
  delay <- 4

  api <- mock_api(
    list(
      human = list(
        ID = level,
        last_boosted_id = last_boosted
      )
    )
  )

  boost_immunity(
    api,
    individuals$human,
    variables$id,
    index,
    level,
    variables$last_boosted_id,
    timestep,
    delay
  )

  mockery::expect_called(api$queue_variable_update, 0)
})


test_that('calculate_treated can handle multiple drugs', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  parameters <- set_clinical_treatment(parameters, .5, c(1, 2), c(.5, .5))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(human = list(infectivity = c(.1, .2, .3))),
    parameters = parameters,
    timestep = 5
  )

  with_mock(
    'malariasimulation:::bernoulli' = mockery::mock(
      seq(3) # Mock seek treatment
    ),
    'malariasimulation:::bernoulli_multi_p' = mockery::mock(
      c(TRUE, TRUE, FALSE) # Mock drug success
    ),
    sample.int = mockery::mock(c(1, 2)),
    calculate_treated(api, individuals$human, states, variables, seq(3))
  )

  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$human,
    states$Tr,
    c(1, 2)
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$infectivity,
    parameters$cd * c(AL_params[[2]], DHC_PQP_params[[2]]),
    c(1, 2)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$drug,
    c(1, 2),
    c(1, 2)
  )

  mockery::expect_args(
    api$queue_variable_update,
    3,
    individuals$human,
    variables$drug_time,
    5,
    c(1, 2)
  )
})
