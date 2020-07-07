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
#' @export
run_simulation <- function(timesteps, parameters = NULL) {
  events <- create_events()
  if (is.null(parameters)) {
    parameters <- get_parameters()
  }
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  mda_events <- create_mda_events(parameters)
  individuals <- create_individuals(
    states,
    variables,
    events,
    parameters,
    mda_events
  )
  create_event_based_processes(individuals, states, variables, events, parameters)
  odes <- parameterise_ode(parameters)
  attach_mda_listeners(individuals$human, states, variables, mda_events)
  individual::simulate(
    individuals = individuals,
    processes = create_processes(
      individuals,
      states,
      variables,
      events,
      parameters,
      odes,
      mda_events
    ),
    end_timestep = timesteps,
    parameters = parameters,
    initialisation = create_setup_process(events)
  )
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
