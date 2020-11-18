initial_immunity <- function(parameter, age) {
  if (length(parameter) == 1) {
    return(rep(parameter, length(age)))
  } else if (length(parameter) == 100) {
    age <- floor(age / 365)
    age[age > 99] <- 99
    return(parameter[age + 1])
  }
  stop('unsupported immunity parameter')
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
#' * Tr - Patients under **Tr**eatment
#'
#' The mosquito states are defined as:
#'
#' * Sm - **S**usceptable **m**osquito
#' * Pm - Extrinsic Incubation **P**eriod for **m**osquitoes
#' * Im - **I**nfectious **m**osquito
#' * Unborn - This is a dummy state to allow for a varying number of mosquitos
#' in our model. Individuals enter the Unborn state when they die and leave when
#' new larvae emerge
#'
#' @param parameters, the model parameters
create_states <- function(parameters) {
  initial_counts <- calculate_initial_counts(parameters)
  states <- list(
    # Human states
    S = individual::State$new("S", initial_counts[[1]]),
    D = individual::State$new("D", initial_counts[[2]]),
    A = individual::State$new("A", initial_counts[[3]]),
    U = individual::State$new("U", initial_counts[[4]]),
    Tr = individual::State$new("Tr", initial_counts[[5]])
  )

  mosquito_counts <- floor(initial_mosquito_counts(parameters, parameters$init_foim))
  n_Unborn <- parameters$mosquito_limit - sum(mosquito_counts[c(4, 5, 6)])

  if (n_Unborn < 0) {
    stop(paste('Mosquito limit not high enough. Short by', -n_Unborn, sep=' '))
  }
  states <- c(
    states,
    # Mosquito states
    Sm      = individual::State$new("Sm", mosquito_counts[[4]]),
    Pm      = individual::State$new("Pm", mosquito_counts[[5]]),
    Im      = individual::State$new("Im", mosquito_counts[[6]]),
    Unborn  = individual::State$new("Unborn", n_Unborn)
  )
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
#' * zeta - Heterogeneity of human individuals
#' * zeta_group - Discretised heterogeneity of human individuals
#' * rtss_vaccinated - The timstep of the last rtss vaccination (-1 if there
#' haven't been any)
#' * rtss_boosted  - The timstep of the last rtss booster (-1 if there
#' haven't been any)
#' * rtss_cs - peak antibodies
#' * rtss_rho - antibody component variable
#' * rtss_ds - short-lived antibody delay variable
#' * rtss_dl - long-lived antibody delay variable
#' * tbv_vaccinated - The timstep of the last tbv vaccination (-1 if there
#' haven't been any
#' * zeta_group - Discretised heterogeneity of human individuals
#' * net_time - The timestep when a net was last put up (-1 if never)
#' * spray_time - The timestep when the house was last sprayed (-1 if never)
#' * infectivity - The onward infectiousness to mosquitos
#' * drug - The last prescribed drug
#' * drug_time - The timestep of the last drug
#'
#' Mosquito variables are: 
#' * variety - The variety of mosquito, this is an ordinal in the range of
#' length(parameters$variety_proportions)
#'
#' @param parameters, model parameters created by `get_parameters`
#' @importFrom stats rexp rnorm
create_variables <- function(parameters) {
  size <- parameters$human_population

  initial_age <- floor(rexp(size, rate=1/parameters$average_age))

  # Define variables
  birth <- individual::Variable$new("birth", -initial_age)
  last_boosted_ib <- individual::Variable$new("last_boosted_ib", rep(-1, size) )
  last_boosted_ica <- individual::Variable$new("last_boosted_ica", rep(-1, size))
  last_boosted_iva <- individual::Variable$new("last_boosted_iva", rep(-1, size))
  last_boosted_id <- individual::Variable$new("last_boosted_id", rep(-1, size))

  is_severe <- individual::Variable$new("is_severe", rep(0, size))

  # Maternal immunity
  icm <- individual::Variable$new(
    "ICM",
    initial_immunity(parameters$init_icm, initial_age)
  )

  ivm <- individual::Variable$new(
    "IVM",
    initial_immunity(parameters$init_ivm, initial_age)
  )

  # Pre-erythoctic immunity
  ib  <- individual::Variable$new(
    "IB",
    initial_immunity(parameters$init_ib, initial_age)
  )
  # Acquired immunity to clinical disease
  ica <- individual::Variable$new(
    "ICA",
    initial_immunity(parameters$init_ica, initial_age)
  )
  # Acquired immunity to severe disease
  iva <- individual::Variable$new(
    "IVA",
    initial_immunity(parameters$init_iva, initial_age)
  )
  # Acquired immunity to detectability
  id <- individual::Variable$new(
    "ID",
    initial_immunity(parameters$init_id, initial_age)
  )

  zeta_norm <- rnorm(size)
  zeta <- individual::Variable$new(
    "zeta",
    exp(
      zeta_norm * sqrt(parameters$sigma_squared) - parameters$sigma_squared/2
    )
  )

  zeta_group <- individual::Variable$new(
    "zeta_group",
    discretise_normal(zeta_norm, parameters$n_heterogeneity_groups)
  )

  # Initialise infectiousness of humans -> mosquitoes
  # NOTE: not yet supporting initialisation of infectiousness of Treated individuals
  infectivity_values <- rep(0, parameters$human_population)
  counts <- calculate_initial_counts(parameters)

  # Calculate the indices of individuals in each infectious state
  diseased <- counts[[1]]:sum(counts[1:2])  # The index of individuals in the D state
  asymptomatic <- sum(counts[1:2]):sum(counts[1:3]) # The index of individuals in the A state
  subpatent <- sum(counts[1:3]):sum(counts[1:4]) # The index of individuals in the U state 

  # Set the initial infectivity values for each individual
  infectivity_values[diseased] <- parameters$cd
  infectivity_values[asymptomatic] <- asymptomatic_infectivity(
    initial_age[asymptomatic],
    initial_immunity(parameters$init_id, initial_age)[asymptomatic],
    parameters
  )
  infectivity_values[subpatent] <- parameters$cu

  # Initialise the infectivity variable
  infectivity <- individual::Variable$new("infectivity", infectivity_values)

  drug <- individual::Variable$new("drug", rep(0, size))
  drug_time <- individual::Variable$new("drug_time", rep(-1, size))

  rtss_vaccinated <- individual::Variable$new("rtss_vaccinated", rep(-1, size))
  rtss_boosted <- individual::Variable$new("rtss_boosted", rep(-1, size))

  rtss_cs <- individual::Variable$new(
    "rtss_cs",
    exp(rnorm(size, parameters$rtss_cs[[1]], parameters$rtss_cs[[2]]))
  )
  rtss_rho <- individual::Variable$new(
    "rtss_rho",
    invlogit(rnorm(size, parameters$rtss_rho[[1]], parameters$rtss_rho[[2]]))
  )
  rtss_ds <- individual::Variable$new(
    "rtss_ds",
    exp(rnorm(size, parameters$rtss_ds[[1]], parameters$rtss_ds[[2]]))
  )

  rtss_dl <- individual::Variable$new(
    "rtss_dl",
    exp(rnorm(size, parameters$rtss_dl[[1]], parameters$rtss_dl[[2]]))
  )

  tbv_vaccinated <- individual::Variable$new("tbv_vaccinated", rep(-1, size))

  # Init vector controls
  net_time <- individual::Variable$new("net_time", rep(-1, size))
  spray_time <- individual::Variable$new("spray_time", rep(-1, size))

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
    zeta = zeta,
    zeta_group = zeta_group,
    infectivity = infectivity,
    drug = drug,
    drug_time = drug_time,
    rtss_vaccinated = rtss_vaccinated,
    rtss_boosted = rtss_boosted,
    rtss_cs = rtss_cs,
    rtss_rho = rtss_rho,
    rtss_ds = rtss_ds,
    rtss_dl = rtss_dl,
    tbv_vaccinated = tbv_vaccinated,
    is_severe = is_severe,
    net_time = net_time,
    spray_time = spray_time
  )

  mosquito_variety <- individual::Variable$new(
    "variety",
    sample(
      seq_along(parameters$variety_proportions),
      parameters$mosquito_limit,
      prob = parameters$variety_proportions,
      replace = TRUE
    )
  )

  variables <- c(
    variables,
    mosquito_variety = mosquito_variety
  )

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
#' @param parameters model parameters
create_individuals <- function(
  states,
  variables,
  events,
  parameters
  ) {
  human <- individual::Individual$new(
    'human',
    states = list(states$S, states$D, states$A, states$U, states$Tr),
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
      variables$zeta,
      variables$zeta_group,
      variables$infectivity,
      variables$drug,
      variables$drug_time,
      variables$rtss_vaccinated,
      variables$rtss_boosted,
      variables$rtss_cs,
      variables$rtss_rho,
      variables$rtss_ds,
      variables$rtss_dl,
      variables$tbv_vaccinated,
      variables$net_time,
      variables$spray_time
    ),
    events = c(
      events$infection,
      events$clinical_infection,
      events$asymptomatic_infection,
      events$subpatent_infection,
      events$recovery,
      events$rtss_vaccination,
      events$rtss_booster,
      events$mda_administer,
      events$smc_enrollment,
      events$smc_administer,
      events$tbv_vaccination,
      events$throw_away_net
    )
  )

  mosquito <- individual::Individual$new(
    'mosquito',
    states=list(
      states$Sm,
      states$Pm,
      states$Im,
      states$Unborn
    ),
    variables = list(variables$mosquito_variety),
    events = list(events$mosquito_infection)
  )

  list(human = human, mosquito = mosquito)
}

calculate_initial_counts <- function(parameters) {
  initial_counts <- round(
    c(
      parameters$s_proportion,
      parameters$d_proportion,
      parameters$a_proportion,
      parameters$u_proportion,
      parameters$t_proportion
    ) * parameters$human_population
  )
  left_over <- parameters$human_population - sum(initial_counts)
  initial_counts[[1]] <- initial_counts[[1]] + left_over
  initial_counts
}
