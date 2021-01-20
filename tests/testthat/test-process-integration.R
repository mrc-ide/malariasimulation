
test_that('create_processes makes valid process functions', {
  parameters <- get_parameters()
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  odes <- parameterise_ode(parameters)
  renderer <- individual::Render$new(1)
  processes <- create_processes(
    renderer,
    variables,
    events,
    parameters,
    odes
  )
  for (process in processes) {
    expect(is.function(process) || inherits(process, 'externalptr'), 'Process is not a function')
  }
})
