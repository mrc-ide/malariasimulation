#' @title Run the simulation
#' @description
#' The main entrypoint for the simulation. run_simulation puts together the
#' model components and runs the malaria simulation. This currently returns a
#' 2D vector (number of humans * timesteps) representing the state of each
#' human at each timestep of the simulation
#'
#' Warning: the return type of this function is likely to change as we figure
#' out what kind of outputs we would like to report from the simulation.
#'
#' @param timesteps, the number of timesteps to run the simulation for
#' @export
run_simulation <- function(timesteps) {
  parameters <- get_parameters()
  states <- create_states()
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables)
  individual::simulate(
    individuals = individuals,
    processes = create_processes(individuals, states, variables, parameters),
    timesteps,
    parameters
  )$render(individuals$human)
}
