#' @title Define model states
#' @description
#' create_states creates the human and mosquito states for the model
#' @param parameters, the model parameters
create_states <- function(parameters) {
  list(
    # Human states
    S       = individual::State$new("S", parameters$human_population),
    I       = individual::State$new("I", 0),
    Treated = individual::State$new("Treated", 0),
    D       = individual::State$new("D", 0),
    A       = individual::State$new("A", 0),
    U       = individual::State$new("U", 0),
    P       = individual::State$new("P", 0),
    # Mosquito states
    E       = individual::State$new("E", 0),
    L       = individual::State$new("L", 0),
    P       = individual::State$new("P", 0),
    Sm      = individual::State$new("Sm", 1),
    Im      = individual::State$new("Im", 0),
    Unborn  = individual::State$new("Unborn", parameters$mosquito_limit - 1)
  )
}

#' @title Define model variables and constants
#' @description
#' create_variables creates the human and mosquito variables and constants for
#' the model. Variables are used to track real data for each individual over
#' time, they are read and updated by processes. Constants remain fixed for each
#' individual over the whole simulation.
#'
#' @param parameters, model parameters created by `get_parameters`
create_variables <- function(parameters) {
  initial_age <- trunc(rexp(parameters$human_population, rate=1/10))

  # Define variables
  age <- individual::Variable$new("age", function(size) initial_age)

  list(
    age = age
  )
}

#' @title Define model individuals
#' @description
#' create_individuals declares the individuals to simulate. It assigns the
#' relevant states and variables to each individual.
#'
#' @param states, available states to assign
#' @param variables, available variables and constants to assign
create_individuals <- function(states, variables) {
  human <- individual::Individual$new(
    'human',
    states=list(states$S, states$I, states$Treated, states$D, states$A, states$U),
    variables = list(variables$age)
  )

  mosquito <- individual::Individual$new(
    'mosquito',
    states=list(states$E, states$L, states$P, states$Sm, states$Im, states$Unborn)
  )

  list(human = human, mosquito = mosquito)
}
