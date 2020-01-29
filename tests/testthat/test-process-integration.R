
test_that('create_processes makes valid process functions', {
  parameters <- get_parameters()
  states <- create_states()
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables)
  processes <- create_processes(individuals, states, variables, parameters)
  for (process in processes) {
    expect(is.function(process), 'Process is not a function')
    expect_equal(length(formals(process)), 3)
  }
})
