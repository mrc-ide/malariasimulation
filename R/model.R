#' @title Run the simulation
#' @description
#' The main entrypoint for the simulation. run_simulation puts together the
#' model components and runs the malaria simulation. This currently returns a
#' dataframe with the number of individuals in each state at each timestep
#'
#' Warning: the columns of the output dataframe is likely to change as we figure
#' out what kind of outputs we would like to report from the simulation.
#'
#' @param timesteps the number of timesteps to run the simulation for
#' @param parameters a named list of parameters to use
#' @param correlations correlation parameters
#' @export
run_simulation <- function(timesteps, parameters = NULL, correlations = NULL) {
  if (is.null(parameters)) {
    parameters <- get_parameters()
  }
  if (is.null(correlations)) {
    correlations <- get_correlation_parameters(parameters)
  }
  variables <- create_variables(parameters)
  events <- create_events(parameters)
  initialise_events(events, variables, parameters)
  renderer <- individual::Render$new(timesteps)
  attach_event_listeners(
    events,
    variables,
    parameters,
    correlations,
    renderer
  )
  vector_models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(vector_models, parameters)
  individual::simulation_loop(
    processes = create_processes(
      renderer,
      variables,
      events,
      parameters,
      vector_models,
      solvers,
      correlations
    ),
    variables = variables,
    events = events,
    timesteps = timesteps
  )
  renderer$to_dataframe()
}

#' @title Run the simulation with repetitions
#'
#' @param timesteps the number of timesteps to run the simulation for
#' @param repetitions n times to run the simulation
#' @param overrides a named list of parameters to use instead of defaults
#' @param parallel execute runs in parallel
#' @export
run_simulation_with_repetitions <- function(
  timesteps,
  repetitions,
  overrides = list(),
  parallel = FALSE
  ) {
  if (parallel) {
    fapply <- parallel::mclapply
  } else {
    fapply <- lapply
  }
  dfs <- fapply(
    seq(repetitions),
    function(repetition) {
      df <- run_simulation(timesteps, overrides)
      df$repetition <- repetition
      df
    }
  )
  do.call("rbind", dfs)
}
