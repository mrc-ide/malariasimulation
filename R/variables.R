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

#' @title Define model variables
#' @description
#' create_variables creates the human and mosquito variables for
#' the model. Variables are used to track real data for each individual over
#' time, they are read and updated by processes
#'
#' The human variables are defined as:
#'
#' * state - the state of a human individual S|D|A|U|Tr
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
#' * mosquito_state - the state of the mosquito, a category Sm|Pm|Im|Unborn
#' * species - the species of mosquito, this is a category gamb|fun|arab
#'
#' @param parameters, model parameters created by `get_parameters`
#' @importFrom stats rexp rnorm
create_variables <- function(parameters) {
  size <- parameters$human_population

  initial_age <- floor(rexp(size, rate=1/parameters$average_age))

  initial_counts <- calculate_initial_counts(parameters)

  # Define variables
  states <- c('S', 'D', 'A', 'U', 'Tr')
  state <- individual::CategoricalVariable$new(
    states,
    rep(states, times = initial_counts)
  )
  birth <- individual::DoubleVariable$new(-initial_age)
  last_boosted_ib <- individual::DoubleVariable$new(rep(-1, size) )
  last_boosted_ica <- individual::DoubleVariable$new(rep(-1, size))
  last_boosted_iva <- individual::DoubleVariable$new(rep(-1, size))
  last_boosted_id <- individual::DoubleVariable$new(rep(-1, size))

  is_severe <- individual::DoubleVariable$new(rep(0, size))

  # Maternal immunity
  icm <- individual::DoubleVariable$new(
    initial_immunity(parameters$init_icm, initial_age)
  )

  ivm <- individual::DoubleVariable$new(
    initial_immunity(parameters$init_ivm, initial_age)
  )

  # Pre-erythoctic immunity
  ib  <- individual::DoubleVariable$new(
    initial_immunity(parameters$init_ib, initial_age)
  )
  # Acquired immunity to clinical disease
  ica <- individual::DoubleVariable$new(
    initial_immunity(parameters$init_ica, initial_age)
  )
  # Acquired immunity to severe disease
  iva <- individual::DoubleVariable$new(
    initial_immunity(parameters$init_iva, initial_age)
  )
  # Acquired immunity to detectability
  id <- individual::DoubleVariable$new(
    initial_immunity(parameters$init_id, initial_age)
  )

  zeta_norm <- rnorm(size)
  zeta <- individual::DoubleVariable$new(
    exp(
      zeta_norm * sqrt(parameters$sigma_squared) - parameters$sigma_squared/2
    )
  )

  zeta_group <- individual::DoubleVariable$new(
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
  infectivity <- individual::DoubleVariable$new(infectivity_values)

  drug <- individual::DoubleVariable$new(rep(0, size))
  drug_time <- individual::DoubleVariable$new(rep(-1, size))

  rtss_vaccinated <- individual::DoubleVariable$new(rep(-1, size))
  rtss_boosted <- individual::DoubleVariable$new(rep(-1, size))

  rtss_cs <- individual::DoubleVariable$new(
    exp(rnorm(size, parameters$rtss_cs[[1]], parameters$rtss_cs[[2]]))
  )
  rtss_rho <- individual::DoubleVariable$new(
    invlogit(rnorm(size, parameters$rtss_rho[[1]], parameters$rtss_rho[[2]]))
  )
  rtss_ds <- individual::DoubleVariable$new(
    exp(rnorm(size, parameters$rtss_ds[[1]], parameters$rtss_ds[[2]]))
  )

  rtss_dl <- individual::DoubleVariable$new(
    exp(rnorm(size, parameters$rtss_dl[[1]], parameters$rtss_dl[[2]]))
  )

  tbv_vaccinated <- individual::DoubleVariable$new(rep(-1, size))

  # Init vector controls
  net_time <- individual::DoubleVariable$new(rep(-1, size))
  spray_time <- individual::DoubleVariable$new(rep(-1, size))

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

  species <- individual::CategoricalVariable$new(
    parameters$species,
    sample(
      parameters$species,
      parameters$mosquito_limit,
      prob = parameters$species_proportions,
      replace = TRUE
    )
  )

  mosquito_counts <- floor(
    initial_mosquito_counts(parameters, parameters$init_foim)
  )

  n_Unborn <- parameters$mosquito_limit - sum(mosquito_counts[c(4, 5, 6)])

  if (n_Unborn < 0) {
    stop(paste('Mosquito limit not high enough. Short by', -n_Unborn, sep=' '))
  }

  mosquito_states <- c('Sm', 'Pm', 'Im', 'Unborn')
  mosquito_state <- individual::CategoricalVariable$new(
    mosquito_states,
    rep(mosquito_states, times = c(mosquito_counts[c(4, 5, 6)], n_Unborn))
  )

  variables <- c(
    variables,
    species = species,
    mosquito_state = mosquito_state
  )

  variables
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
