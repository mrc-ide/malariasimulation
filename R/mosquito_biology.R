#' @title Calculate equilibrium solution for vector counts
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param species the index of the species to find the equilibrium for
#' @param foim equilibrium foim
#' @param m the total number of female adult mosquitos
#' @noRd
initial_mosquito_counts <- function(parameters, species, foim, m) {
  omega <- calculate_omega(parameters, species)
  mum <- parameters$mum[[species]]
  n_E <- 2 * omega * mum * parameters$dl * (
    1. + parameters$dpl * parameters$mup
  ) * m

  n_L <- 2 * mum * parameters$dl * (
    1. + parameters$dpl * parameters$mup
  ) * m

  n_P <- 2 * parameters$dpl * mum * m

  n_Sm <- m * mum / (foim + mum)

  incubation_survival <- exp(-mum * parameters$dem)

  n_Pm <- m * foim / (foim + mum) * (
    1. - incubation_survival
  )

  n_Im <- m * foim / (foim + mum) * incubation_survival

  c(n_E, n_L, n_P, n_Sm, n_Pm, n_Im)
}

#' @title Calculate omega value
#' @description useful value for calculating equilibrium solutions for vectors
#' taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param species the index of the species to calculate for
#' @noRd
calculate_omega <- function(parameters, species) {
  sub_omega <- parameters$gamma * parameters$ml / parameters$me - (
    parameters$del / parameters$dl
  ) + (
    (parameters$gamma - 1) * parameters$ml * parameters$del
  )

  mum <- parameters$mum[[species]]

  beta <- eggs_laid(
    parameters$beta,
    mum,
    parameters$blood_meal_rates[[species]]
  )

  -.5 * sub_omega + sqrt(
    .25 * sub_omega**2 +
      .5 * parameters$gamma * beta * parameters$ml * parameters$del /
      (parameters$me * mum * parameters$dl * (
        1. + parameters$dpl * parameters$mup
      ))
  )
}

#' @title Calculate the vector carrying capacity
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param m number of adult mosquitoes
#' @param species index of the species to calculate for
calculate_carrying_capacity <- function(parameters, m, species) {
  omega <- calculate_omega(parameters, species)

  m * 2 * parameters$dl * parameters$mum[[species]] * (
    1. + parameters$dpl * parameters$mup
  ) * parameters$gamma * (omega + 1) / (
    omega / (parameters$ml * parameters$del) - (
      1. / (parameters$ml * parameters$dl)
    ) - 1.
  )
}

#' @title Calculate the mean rainfall throughout the year
#' @param parameters model parameters
#' @noRd
calculate_R_bar <- function(parameters) {
  mean(vnapply(1:365, function(t) rainfall(
		t,
    parameters$g0,
    parameters$g,
    parameters$h,
    parameters$rainfall_floor
	)))
}

#' @title Calculate equilibrium total_M from parameters
#'
#' @param parameters to work from
#' @param EIR equilibrium to use, bites per person per year
#' @importFrom stats weighted.mean
#' @noRd
equilibrium_total_M <- function(parameters, EIR) {
  if (EIR == 0) {
    return(0)
  }
  if (parameters$init_foim == 0) {
    stop('init_foim must be > 0 to calculate a non-zero equilibrium total_M')
  }
  mum <- weighted.mean(parameters$mum, parameters$species_proportions)
  total_daily_eir <- EIR * parameters$human_population / 365
  lifetime <- parameters$init_foim * exp(-mum * parameters$dem) / (
    parameters$init_foim + mum
  )
  total_daily_eir / sum(
    parameters$species_proportions * parameters$blood_meal_rates * parameters$Q0 * lifetime
  )
}

#' @title Calculate the yearly offset (in timesteps) for the peak mosquito
#' season
#'
#' @param parameters to work from
#' @export
peak_season_offset <- function(parameters) {
  if (!parameters$model_seasonality) {
    return(0)
  }
  which.max(vnapply(seq(365), function(t) {
    rainfall(
      t,
      parameters$g0,
      parameters$g,
      parameters$h,
      parameters$rainfall_floor
    )
  }))[[1]]
}

#' @title Calculate the death rate of mosquitoes given interventions
#'
#' @param f the feeding rate for this species of mosquito
#' @param W the mean probability that a mosquito feeds and survives
#' @param Z the mean probability that a mosquito is repelled
#' @param Z the mean probability that a mosquito is repelled
#' @noRd
death_rate <- function(f, W, Z, species, parameters) {
  mum <- parameters$mum[[species]]
  p1_0 <- exp(-mum * parameters$foraging_time[[species]])
  gonotrophic_cycle <- get_gonotrophic_cycle(species, parameters)
  p2 <- exp(-mum * gonotrophic_cycle)
  p1 <- p1_0 * W / (1 - Z * p1_0)
  -f * log(p1 * p2)
}

get_gonotrophic_cycle <- function(v, parameters) {
  f <- parameters$blood_meal_rates[[v]]
  gonotrophic_cycle <- 1 / f - parameters$foraging_time[[v]]
}

#' @title Update the individual mosquito model after biting
#'
#' @param variables a list of variables in this simulation
#' @param foim force of infection towards mosquitoes
#' @param events events in the simulation
#' @param species the index of the species to calculate for
#' @param susceptible_species the indices of susceptible mosquitos of the
#' species
#' @param adult_species the indices of adult mosquitos of the species
#' @param mu the death rate of the current species
#' @param parameters the model parameters
#' @param timestep the current timestep
#' @noRd
biting_effects_individual <- function(
    variables,
    foim,
    events,
    species,
    susceptible_species,
    adult_species,
    mu,
    parameters,
    timestep
  ) {
  # deal with mosquito infections
  target <- sample_bitset(susceptible_species, foim)
  variables$mosquito_state$queue_update('Pm', target)
  events$mosquito_infection$schedule(
    target,
    log_uniform(target$size(), parameters$dem)
  )

  # deal with mosquito deaths
  died <- sample_bitset(adult_species, mu)

  events$mosquito_death$schedule(died, 0)
}
#' @title Mosquito emergence process
#' @description Move mosquitos from NonExistent to Sm in line with the number of
#' pupals in the ODE models
#'
#' @param solvers a list of solver objects for each species of mosquito
#' @param state the variable for the mosquito state
#' @param species the variable for the mosquito species
#' @param species_names a character vector of species names for each solver
#' @param dpl the delay for pupal growth (in timesteps)
#' @noRd
create_mosquito_emergence_process <- function(
  solvers,
  state,
  species,
  species_names,
  dpl
  ) {
  rate <- .5 * 1 / dpl
  function(timestep) {
    p_counts <- vnapply(
      solvers,
      function(solver) {
        solver_get_states(solver)[[ODE_INDICES[['P']]]]
      }
    )
    n <- sum(p_counts) * rate
    available <- state$get_size_of('NonExistent')
    if (n > available) {
      stop(paste0(
        'Not enough mosquitoes (short by ',
        n - available,
        '). Please raise parameters$mosquito_limit. ',
        'If you have used parameterise_mosquito_equilibrium,',
        'your seasonality parameters lead to more mosquitoes than expected.'
      ))
    }
    non_existent <- state$get_index_of('NonExistent')
    latest <- 1
    for (i in seq_along(species_names)) {
      to_hatch <- p_counts[[i]] * rate
      hatched <- bitset_at(non_existent, seq(latest, latest + to_hatch))
      state$queue_update('Sm', hatched)
      species$queue_update(species_names[[i]], hatched)
      latest <- latest + to_hatch + 1
    }
  }
}
