#' @title Define model states
#' @description
#' create_states creates the human and mosquito states for the model
#' 
#' The human states are defined as:
#' 
#' * S - **S**usceptable to infection
#' * I - Liver-stage **I**nfection, these individuals are waiting to develop to A
#' or D levels of infection
#' * D - **D**isease individuals exhibit "clinical" or "severe" disease
#' * A - **A**symptomatic individuals no longer exhibit symptoms
#' * U - S**u**bpatent infectious patients are still infectious to mosquitos
#'
#' The mosquito states are defined as:
#'
#' * E - **E**arly larval stage
#' * L - **L**ate larval stage
#' * P - **P**upal
#' * Sm - **S**usceptable **m**osquito
#' * Im - **I**nfectious **m**osquito
#' * Unborn - This is a dummy state to allow for a varying number of mosquitos
#' in our model. Individuals enter the Unborn state when they die and leave when
#' new larvae emerge
#'
#' @param parameters, the model parameters
create_states <- function(parameters) {
  list(
    # Human states
    S       = individual::State$new("S", parameters$human_population),
    I       = individual::State$new("I", 0),
    D       = individual::State$new("D", 0),
    A       = individual::State$new("A", 0),
    U       = individual::State$new("U", 0),
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
#' The human variables are defined as:
#'
#' * age - an integer representing the number of years this individual has been
#' alive
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
    states=list(states$S, states$I, states$D, states$A, states$U),
    variables = list(variables$age)
  )

  mosquito <- individual::Individual$new(
    'mosquito',
    states=list(states$E, states$L, states$P, states$Sm, states$Im, states$Unborn)
  )

  list(human = human, mosquito = mosquito)
}
