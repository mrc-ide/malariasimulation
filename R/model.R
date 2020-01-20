#' @description
#' 
#' Run malaria simulation
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
