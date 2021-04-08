
test_that('create_processes makes valid process functions', {
  parameters <- get_parameters()
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  vector_models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(vector_models, parameters)
  renderer <- individual::Render$new(1)
  processes <- create_processes(
    renderer,
    variables,
    events,
    parameters,
    vector_models,
    solvers
  )
  for (process in processes) {
    expect(is.function(process) || inherits(process, 'externalptr'), 'Process is not a function')
  }
})

test_that('attach_event_listeners makes valid listeners', {
  parameters <- get_parameters()
  events <- create_events(parameters)
  correlations <- get_correlation_parameters(parameters)
  variables <- create_variables(parameters)
  renderer <- individual::Render$new(2)
  attach_event_listeners(events, variables, parameters, correlations, renderer)

  for (event in events) {
    for (listener in event$.listeners) {
      expect(
        (is.function(listener) && length(args(listener)) >= 1 && length(args(listener)) <= 2)  
          || inherits(process, 'externalptr'),
        'Listener is not valid'
      )
    }
  }
})
