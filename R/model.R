#' @title Run the simulation
#' @description
#' The main entrypoint for the simulation. run_simulation puts together the
#' model components and runs the malaria simulation. This currently returns a
#' dataframe with the number of individuals in each state at each timestep
#'
#' Warning: the return type of this function is likely to change as we figure
#' out what kind of outputs we would like to report from the simulation.
#'
#' @param timesteps, the number of timesteps to run the simulation for
#' @export
run_simulation <- function(timesteps, overrides = list()) {
  parameters <- get_parameters(overrides)
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables)
  individual::simulate(
    individuals = individuals,
    processes = create_processes(individuals, states, variables, parameters),
    end_timestep = timesteps,
    parameters = parameters,
    custom_renderers = create_renderers(individuals, states, variables, parameters)
  )
}
