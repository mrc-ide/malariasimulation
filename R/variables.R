#' @title Define model variables
#' @description
#' create_variables creates the human and mosquito variables for
#' the model. Variables are used to track real data for each individual over
#' time, they are read and updated by processes
#'
#' The human variables are defined as:
#'
#' * state - the state of a human individual S|D|A|U|Tr|NonExistent
#' * birth - an integer representing the timestep when this individual was born
#' * last_boosted_* - the last timestep at which this individual's immunity was
#' boosted for tracking grace periods in the boost of immunity
#' * is_severe - yes or no if the individual currently has severe malaria
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
#' * net_time - The timestep when a net was last put up (-1 if never)
#' * spray_time - The timestep when the house was last sprayed (-1 if never)
#' * infectivity - The onward infectiousness to mosquitos
#' * drug - The last prescribed drug
#' * drug_time - The timestep of the last drug
#'
#' Mosquito variables are: 
#' * mosquito_state - the state of the mosquito, a category Sm|Pm|Im|NonExistent
#' * species - the species of mosquito, this is a category gamb|fun|arab
#'
#' @param parameters, model parameters created by `get_parameters`
#' @noRd
#' @importFrom stats rexp rnorm
create_variables <- function(parameters) {
  size <- parameters$human_population

  initial_age <- calculate_initial_ages(parameters)

  if (parameters$enable_heterogeneity) {
    quads <- statmod::gauss.quad.prob(
      parameters$n_heterogeneity_groups,
      dist='normal'
    )
    groups <- sample.int(
      parameters$n_heterogeneity_groups,
      size,
      replace = TRUE,
      prob = quads$weights
    )
    zeta_norm <- quads$nodes[groups]
    zeta <- individual::DoubleVariable$new(
      calculate_zeta(zeta_norm, parameters)
    )
    zeta_group <- individual::CategoricalVariable$new(
      to_char_vector(seq(parameters$n_heterogeneity_groups)),
      to_char_vector(groups)
    )
    if (!is.null(parameters$init_EIR)) {
      eq <- calculate_eq(quads$nodes, parameters)
    } else {
      eq <- NULL
    }
  } else {
    zeta <- individual::DoubleVariable$new(rep(1, size))
    groups <- rep(1, size)
    zeta_group <- individual::CategoricalVariable$new(
      to_char_vector(seq(parameters$n_heterogeneity_groups)),
      to_char_vector(groups)
    )
    if (!is.null(parameters$init_EIR)) {
      eq <-	list(
        malariaEquilibrium::human_equilibrium_no_het(
          parameters$init_EIR,
          sum(get_treatment_coverages(parameters, 0)),
          parameters$eq_params,
          EQUILIBRIUM_AGES
        )
      )
    } else {
      eq <- NULL
    }
  }

  states <- c('S', 'D', 'A', 'U', 'Tr', 'NonExistent')
  state <- individual::CategoricalVariable$new(
    states,
    initial_state(parameters, initial_age, groups, eq)
  )
  birth <- individual::IntegerVariable$new(-initial_age)
  last_boosted_ib <- individual::DoubleVariable$new(rep(-1, size))
  last_boosted_ica <- individual::DoubleVariable$new(rep(-1, size))
  last_boosted_iva <- individual::DoubleVariable$new(rep(-1, size))
  last_boosted_id <- individual::DoubleVariable$new(rep(-1, size))

  is_severe <- individual::CategoricalVariable$new(
    c('yes', 'no'),
    rep('no', size)
  )

  # Maternal immunity
  icm <- individual::DoubleVariable$new(
    initial_immunity(
      parameters$init_icm,
      initial_age,
      groups,
      eq,
      parameters,
      'ICM'
    )
  )

  ivm <- individual::DoubleVariable$new(
    initial_immunity(parameters$init_ivm, initial_age)
  )

  # Pre-erythoctic immunity
  ib  <- individual::DoubleVariable$new(
    initial_immunity(
      parameters$init_ib,
      initial_age,
      groups,
      eq,
      parameters,
      'IB'
    )
  )
  # Acquired immunity to clinical disease
  ica <- individual::DoubleVariable$new(
    initial_immunity(
      parameters$init_ica,
      initial_age,
      groups,
      eq,
      parameters,
      'ICA'
    )
  )
  # Acquired immunity to severe disease
  iva <- individual::DoubleVariable$new(
    initial_immunity(parameters$init_iva, initial_age)
  )
  # Acquired immunity to detectability
  id <- individual::DoubleVariable$new(
    initial_immunity(
      parameters$init_id,
      initial_age,
      groups,
      eq,
      parameters,
      'ID'
    )
  )

  # Initialise infectiousness of humans -> mosquitoes
  # NOTE: not yet supporting initialisation of infectiousness of Treated individuals
  infectivity_values <- rep(0, parameters$human_population)
  counts <- calculate_initial_counts(parameters)

  # Calculate the indices of individuals in each infectious state
  diseased <- state$get_index_of('D')$to_vector()
  asymptomatic <- state$get_index_of('A')$to_vector()
  subpatent <- state$get_index_of('U')$to_vector()

  # Set the initial infectivity values for each individual
  infectivity_values[diseased] <- parameters$cd
  infectivity_values[asymptomatic] <- asymptomatic_infectivity(
    initial_age[asymptomatic],
    id$get_values(asymptomatic),
    parameters
  )
  infectivity_values[subpatent] <- parameters$cu

  # Initialise the infectivity variable
  infectivity <- individual::DoubleVariable$new(infectivity_values)

  drug <- individual::IntegerVariable$new(rep(0, size))
  drug_time <- individual::IntegerVariable$new(rep(-1, size))

  rtss_vaccinated <- individual::IntegerVariable$new(rep(-1, size))
  rtss_boosted <- individual::IntegerVariable$new(rep(-1, size))

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
  net_time <- individual::IntegerVariable$new(rep(-1, size))
  spray_time <- individual::IntegerVariable$new(rep(-1, size))

  variables <- list(
    state = state,
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

  # Add variables for individual mosquitoes
  if (parameters$individual_mosquitoes) {
    species_values <- NULL
    state_values <- NULL
    n_initialised <- 0
    for (i in seq_along(parameters$species)) {
      mosquito_counts <- floor(
        initial_mosquito_counts(
          parameters,
          i,
          parameters$init_foim,
          parameters$total_M * parameters$species_proportions[[i]]
        )
      )

      species_M <- sum(mosquito_counts[ADULT_ODE_INDICES])

      if (species_M > 0) {
        if (length(species_values) > parameters$mosquito_limit) {
          stop('Mosquito limit not high enough')
        }

        species_values <- c(
          species_values,
          rep(parameters$species[[i]], species_M)
        )
        state_values <- c(
          state_values,
          rep(
            c('Sm', 'Pm', 'Im'),
            times = mosquito_counts[ADULT_ODE_INDICES]
          )
        )
      }
    }

    # fill excess mosquitoes
    excess <- parameters$mosquito_limit - length(species_values)
    species_values <- c(
      species_values,
      rep(parameters$species[[1]], excess)
    )
    state_values <- c(state_values, rep('NonExistent', excess))

    # initialise variables
    species <- individual::CategoricalVariable$new(
      parameters$species,
      species_values
    )
    mosquito_state <- individual::CategoricalVariable$new(
      c('Sm', 'Pm', 'Im', 'NonExistent'),
      state_values
    )
    variables <- c(
      variables,
      species = species,
      mosquito_state = mosquito_state
    )
  }

  variables
}

# =========
# Utilities
# =========

initial_immunity <- function(
  parameter,
  age,
  groups = NULL,
  eq = NULL,
  parameters = NULL,
  eq_name = NULL
  ) {
  if (!is.null(eq)) {
    age <- age / 365
    return(vnapply(
      seq_along(age),
      function(i) {
        g <- groups[[i]]
        a <- age[[i]]
        eq[[g]][which.max(a < eq[[g]][, 'age']), eq_name]
      }
    ))
  }
  rep(parameter, length(age))
}

initial_state <- function(parameters, age, groups, eq) {
  ibm_states <- c('S', 'A', 'D', 'U', 'Tr')
  if (!is.null(eq)) {
    eq_states <- c('S', 'A', 'D', 'U', 'T')
    age <- age / 365
    return(vcapply(
      seq_along(age),
      function(i) {
        g <- groups[[i]]
        a <- age[[i]]
        sample(
          ibm_states,
          size = 1,
          prob = eq[[g]][which.max(a < eq[[g]][, 'age']), eq_states]
        )
      }
    ))
  }
  rep(ibm_states, times = calculate_initial_counts(parameters))
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

calculate_eq <- function(het_nodes, parameters) {
  ft <- sum(get_treatment_coverages(parameters, 0))
	lapply(
		het_nodes,
		function(n) {
			malariaEquilibrium::human_equilibrium_no_het(
				parameters$init_EIR * calculate_zeta(n, parameters),
				ft,
				parameters$eq_params,
        EQUILIBRIUM_AGES
			)
		}
	)
}

calculate_zeta <- function(zeta_norm, parameters) {
  exp(
    zeta_norm * sqrt(parameters$sigma_squared) - parameters$sigma_squared/2
  )
}

calculate_initial_ages <- function(parameters) {
  # check if we've set up a custom demography
  if (is.null(parameters$demography_agegroups)) {
    return(round(rexp(
      parameters$human_population,
      rate=1/parameters$average_age
    )))
  }

  age_high <- parameters$demography_agegroups
  n_age <- length(age_high)
  birthrate <- parameters$demography_birthrates[[1]]
  deathrate <- parameters$demography_deathrates[,1]

  # r[i] can be thought of as the rate of ageing in this age group, i.e.
  # 1/r[i] is the duration of this group
  age_width <- diff(c(0, age_high))
  aging_rate <- 1 / age_width
  aging_rate[[n_age]] <- 0

  prop <- rep(0, n_age)
  for (i in seq_along(parameters$demography_agegroups)) {
    # calculate proportions in each age group
    # birth rate inflow / (aging out + deathrate)
    if (i == 1) {
      prop[i] <- birthrate / (aging_rate[[i]] + deathrate[[i]])
    } else {
      prop[i] <- prop[[i-1]] * aging_rate[[i-1]] / (aging_rate[i] + deathrate[[i]])
    }
  }

  sampled_age_high <- sample(
    age_high,
    parameters$human_population,
    replace = TRUE,
    prob = prop
  )

  group_dist <- runif(0, 1, parameters$human_population)

  sampled_age_high - group_dist * age_width
}
