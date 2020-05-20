
test_that('create_processes makes valid process functions', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables)
  processes <- create_processes(individuals, states, variables, events, parameters)
  for (process in processes) {
    expect(is.function(process) || inherits(process, 'externalptr'), 'Process is not a function')
  }
})
