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
#' * last_bitten - the last timestep at which this individual was bitten, used
#' for tracking grace periods in the boost of immunity
#' * last_infected - the last timestep at which this individual was infected, used
#' for tracking grace periods in the boost of immunity
#' * infection_schedule - the last timestep at which this individual was infected, used
#' for scheduling state transitions to I
#' * asymptomatic_infection_schedule - the last timestep at which this
#' individual was infected, used for scheduling state transitions to A
#' * is_severe - a binary indicator (0 or 1) for if the individual currently has
#' severe malaria
#' * ICM - Maternal immunity to clinical disease
#' * IVM - Maternal immunity to severe disease
#' * IB  - Pre-erythoctic immunity
#' * ICA  - Acquired immunity to clinical disease
#' * IVA  - Acquired immunity to severe disease
#' * ID - Acquired immunity to detectability
#' * xi - Heterogeneity of human individuals
#' * variety - The variety of mosquito, either 1, 2 or 3. These are related to
#' blood meal rate parameters av1, av2 and av3
#'
#' @param parameters, model parameters created by `get_parameters`
create_variables <- function(parameters) {
  initial_age <- trunc(rexp(parameters$human_population, rate=1/10))

  # Define variables
  age <- individual::Variable$new("age", function(size) initial_age)
  last_bitten <- individual::Variable$new("last_bitten", function(size) { rep(-1, size) })
  last_infected <- individual::Variable$new("last_infected", function(size) { rep(-1, size) })
  infection_schedule <- individual::Variable$new(
    "infection_schedule",
    function(size) { rep(-1, size) }
  )
  asymptomatic_infection_schedule <- individual::Variable$new(
    "asymptomatic_infection_schedule",
    function(size) { rep(-1, size) }
  )
  is_severe <- individual::Variable$new(
    "is_severe",
    function(size) { rep(0, size) }
  )

  # Maternal immunity
  icm <- individual::Variable$new(
    "ICM",
    function(size) {
      first_immunity <- 1
      t <- initial_age * 365 / parameters$days_per_timestep
      first_immunity * exp(-(t * parameters$rm))
    }
  )

  ivm <- individual::Variable$new(
    "IVM",
    function(size) {
      first_immunity <- 1
      t <- initial_age * 365 / parameters$days_per_timestep
      first_immunity * exp(-(t * parameters$rm))
    }
  )

  # Pre-erythoctic immunity
  ib  <- individual::Variable$new("IB", function(size) { rep(0, size) })
  # Acquired immunity to clinical disease
  ica <- individual::Variable$new("ICA", function(size) { rep(0, size) })
  # Acquired immunity to severe disease
  iva <- individual::Variable$new("IVA", function(size) { rep(0, size) })
  # Acquired immunity to detectability
  id <- individual::Variable$new("ID", function(size) { rep(0, size) })

  xi_values <- rlnorm(parameters$human_population, -parameters$sigma_squared/2,parameters$sigma_squared)
  xi <- individual::Constant$new(
    "xi",
    function(n) {
      xi_values
    }
  )

  xi_group <- individual::Constant$new(
    "xi_group",
    function(n) {
      discretise(xi_values, n_heterogeneity_groups)
    }
  )

  mosquito_variety <- individual::Constant$new(
    "variety",
    function(n) {
      p <- runif(n)
      v <- rep(0, n)
      v[which(p > .5)] <- 1
      v[which(p > .2 & p < .5)] <- 2
      v[which(p < .2)] <- 3
      v
    }
  )
  list(
    age = age,
    last_bitten = last_bitten,
    last_infected = last_infected,
    icm = icm,
    ivm = ivm,
    ib = ib,
    ica = ica,
    iva = iva,
    id = id,
    xi = xi,
    xi_group = xi_group,
    mosquito_variety = mosquito_variety,
    infection_schedule = infection_schedule,
    asymptomatic_infection_schedule = asymptomatic_infection_schedule,
    is_severe = is_severe
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
    variables = list(
      variables$age,
      variables$last_bitten,
      variables$last_infected,
      variables$ib,
      variables$ica,
      variables$iva,
      variables$id,
      variables$icm,
      variables$ivm,
      variables$infection_schedule,
      variables$asymptomatic_infection_schedule,
      variables$is_severe
    ),
    constants = list(variables$xi)
  )

  mosquito <- individual::Individual$new(
    'mosquito',
    states=list(states$E, states$L, states$P, states$Sm, states$Im, states$Unborn),
    constants = list(variables$mosquito_variety)
  )

  list(human = human, mosquito = mosquito)
}
