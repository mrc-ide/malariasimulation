
test_that('create_processes makes valid process functions', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  odes <- parameterise_ode(parameters)
  processes <- create_processes(
    individuals,
    states,
    variables,
    events,
    parameters,
    odes
  )
  for (process in processes) {
    expect(is.function(process) || inherits(process, 'externalptr'), 'Process is not a function')
  }
})

test_that('Infection listener updates infectivity correctly', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  create_event_based_processes(individuals, states, variables, events, parameters)

  api <- mock_api(list(), parameters = parameters)

  events$infection$listeners[[2]](api, c(2, 4))
  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$infectivity,
    parameters$cd,
    c(2, 4)
  )
})

test_that('Infection listener doesn\'t call for empty targets', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  create_event_based_processes(individuals, states, variables, events, parameters)

  api <- mock_api(list(), parameters = parameters)

  events$infection$listeners[[2]](api, numeric(0))
  mockery::expect_called(api$queue_variable_update, 0)
})

test_that('Asymptomatic infection listener updates infectivity correctly', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  create_event_based_processes(individuals, states, variables, events, parameters)

  api <- mock_api(list(), parameters = parameters)

  with_mock(
    'malariasimulation:::asymptomatic_infectivity' = mockery::mock(c(.2, .5)),
    events$asymptomatic_infection$listeners[[2]](api, c(2, 4))
  )
  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$infectivity,
    c(0.2, 0.5),
    c(2, 4)
  )
})
