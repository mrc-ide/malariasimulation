initial_immunity <- function(parameter, age) {
  if (length(parameter) == 1) {
    return(rep(parameter, length(age)))
  } else if (length(parameter) == 100) {
    age <- trunc(age / 365)
    age[age > 99] <- 99
    return(parameter[age + 1])
  }
  stop('unsupported param')
}

#' @title Define model states
#' @description
#' create_states creates the human and mosquito states for the model
#' 
#' The human states are defined as:
#' 
#' * S - **S**usceptable to infection
#' progressing to the T or D state
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
#' * Pm - Extrinsic Incubation **P**eriod for **m**osquitoes
#' * Im - **I**nfectious **m**osquito
#' * Unborn - This is a dummy state to allow for a varying number of mosquitos
#' in our model. Individuals enter the Unborn state when they die and leave when
#' new larvae emerge
#'
#' @param parameters, the model parameters
create_states <- function(parameters) {
  initial_counts <- vnapply(
    c(
      parameters$s_proportion,
      parameters$d_proportion,
      parameters$a_proportion,
      parameters$u_proportion
    ),
    function(p) round(parameters$human_population * p)
  )
  left_over <- parameters$human_population - sum(initial_counts)
  initial_counts[[1]] <- initial_counts[[1]] + left_over

  states <- list(
    # Human states
    S = individual::State$new(
      "S",
      initial_counts[[1]]
    ),
    D = individual::State$new(
      "D",
      initial_counts[[2]]
    ),
    A = individual::State$new(
      "A",
      initial_counts[[3]]
    ),
    U = individual::State$new(
      "U",
      initial_counts[[4]]
    )
  )

  if (!parameters$vector_ode) {
    mosquito_counts <- initial_mosquito_counts(parameters, parameters$init_foim)
    n_Unborn <- parameters$mosquito_limit - sum(mosquito_counts)

    if (n_Unborn < 0) {
      stop(paste('Mosquito limit not high enough. Short', n_Unborn, sep=' '))
    }
    states <- c(
      states,
      # Mosquito states
      E       = individual::State$new("E", mosquito_counts[[1]]),
      L       = individual::State$new("L", mosquito_counts[[2]]),
      P       = individual::State$new("P", mosquito_counts[[3]]),
      Sm      = individual::State$new("Sm", mosquito_counts[[4]]),
      Pm      = individual::State$new("Pm", mosquito_counts[[5]]),
      Im      = individual::State$new("Im", mosquito_counts[[6]]),
      Unborn  = individual::State$new("Unborn", n_Unborn)
    )
  }
  states
}

#' @title Define model variables
#' @description
#' create_variables creates the human and mosquito variables for
#' the model. Variables are used to track real data for each individual over
#' time, they are read and updated by processes
#'
#' The human variables are defined as:
#'
#' * birth - an integer representing the timestep when this individual was born
#' * last_boosted_* - the last timestep at which this individual's immunity was
#' boosted for tracking grace periods in the boost of immunity
#' * is_severe - a binary indicator (0 or 1) for if the individual currently has
#' severe malaria
#' * ICM - Maternal immunity to clinical disease
#' * IVM - Maternal immunity to severe disease
#' * IB  - Pre-erythoctic immunity
#' * ICA  - Acquired immunity to clinical disease
#' * IVA  - Acquired immunity to severe disease
#' * ID - Acquired immunity to detectability
#' * xi - Heterogeneity of human individuals
#' * xi_group - Discretised heterogeneity of human individuals
#' * variety - The variety of mosquito, either 1, 2 or 3. These are related to
#' blood meal rate parameter
#'
#' @param parameters, model parameters created by `get_parameters`
create_variables <- function(parameters) {
  initial_age <- trunc(
    rexp(
         parameters$human_population,
         rate=1/parameters$average_age
    )
  )

  # Define variables
  birth <- individual::Variable$new("birth", function(size) -initial_age)
  last_boosted_ib <- individual::Variable$new("last_boosted_ib", function(size) { rep(-1, size) })
  last_boosted_ica <- individual::Variable$new("last_boosted_ica", function(size) { rep(-1, size) })
  last_boosted_iva <- individual::Variable$new("last_boosted_iva", function(size) { rep(-1, size) })
  last_boosted_id <- individual::Variable$new("last_boosted_id", function(size) { rep(-1, size) })
  is_severe <- individual::Variable$new(
    "is_severe",
    function(size) { rep(0, size) }
  )

  # Maternal immunity
  icm <- individual::Variable$new(
    "ICM",
    function(size) {
      first_immunity <- parameters$init_icm
      t <- initial_age * 365 / parameters$days_per_timestep
      first_immunity * exp(-(t * parameters$rm))
    }
  )

  ivm <- individual::Variable$new(
    "IVM",
    function(size) {
      first_immunity <- parameters$init_ivm
      t <- initial_age * 365 / parameters$days_per_timestep
      first_immunity * exp(-(t * parameters$rm))
    }
  )

  # Pre-erythoctic immunity
  ib  <- individual::Variable$new(
    "IB",
    function(size) initial_immunity(parameters$init_ib, initial_age)
  )
  # Acquired immunity to clinical disease
  ica <- individual::Variable$new(
    "ICA",
    function(size) initial_immunity(parameters$init_ica, initial_age)
  )
  # Acquired immunity to severe disease
  iva <- individual::Variable$new(
    "IVA",
    function(size) initial_immunity(parameters$init_iva, initial_age)
  )
  # Acquired immunity to detectability
  id <- individual::Variable$new(
     "ID",
    function(size) initial_immunity(parameters$init_id, initial_age)
  )

  xi_values <- exp(rnorm(
    parameters$human_population,
    -parameters$sigma_squared/2,
    sqrt(parameters$sigma_squared)
  ))
  xi <- individual::Variable$new(
    "xi",
    function(n) {
      xi_values
    }
  )

  xi_group <- individual::Variable$new(
    "xi_group",
    function(n) {
      discretise(xi_values, parameters$n_heterogeneity_groups)
    }
  )

  variables <- list(
    birth = birth,
    last_boosted_ib = last_boosted_ib,
    last_boosted_ica = last_boosted_ica,
    last_boosted_iva = last_boosted_iva,
    last_boosted_id = last_boosted_id,
    icm = icm,
    ivm = ivm,
    ib = ib,
    ica = ica,
    iva = iva,
    id = id,
    xi = xi,
    xi_group = xi_group,
    is_severe = is_severe
  )

  if (!parameters$vector_ode) {
    mosquito_variety <- individual::Variable$new(
      "variety",
      function(n) {
        p <- runif(n)
        v <- rep(0, n)
        bottom <- 0
        for (i in seq_along(parameters$variety_proportions)) {
          if (i == 1) {
            match_bottom <- (p >= bottom)
          } else {
            match_bottom <- (p > bottom)
          }
          if (i == length(parameters$variety_proportions)) {
            match_top <- (p <= bottom + parameters$variety_proportions[[i]])
          } else {
            match_top <- (p < bottom + parameters$variety_proportions[[i]])
          }
          v[match_bottom & match_top] <- i
          bottom <- bottom + parameters$variety_proportions[[i]]
        }
        v
      }
    )

    variables <- c(
      variables,
      mosquito_variety = mosquito_variety
    )
  }
  variables
}

#' @title Define model individuals
#' @description
#' create_individuals declares the individuals to simulate. It assigns the
#' relevant states and variables to each individual.
#'
#' @param states available states to assign
#' @param variables available variables to assign
#' @param events available events to assign
create_individuals <- function(states, variables, events, parameters) {
  human <- individual::Individual$new(
    'human',
    states = list(states$S, states$D, states$A, states$U),
    variables = list(
      variables$birth,
      variables$last_boosted_ib,
      variables$last_boosted_ica,
      variables$last_boosted_iva,
      variables$last_boosted_id,
      variables$ib,
      variables$ica,
      variables$iva,
      variables$id,
      variables$icm,
      variables$ivm,
      variables$is_severe,
      variables$xi,
      variables$xi_group
    ),
    events = list(
      events$infection,
      events$asymptomatic_infection
    )
  )

  if (parameters$vector_ode) {
    return(list(human = human))
  }

  mosquito <- individual::Individual$new(
    'mosquito',
    states=list(
      states$E,
      states$L,
      states$P,
      states$Sm,
      states$Pm,
      states$Im,
      states$Unborn
    ),
    variables = list(variables$mosquito_variety),
    events = list(
      events$mosquito_infection
    )
  )

  list(human = human, mosquito = mosquito)
}
