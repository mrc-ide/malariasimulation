create_transmission_mixer <- function(
  variables,
  parameters,
  lagged_eir,
  lagged_infectivity,
  mixing_tt,
  mixing,
  p_captured_tt,
  p_captured,
  p_success
  ) {
  function (timestep) {
    n_pops <- length(variables)
    # calculate all the mixing values once per timestep
    rdt_negative <- vnapply(
      seq(n_pops),
      function(i) {
        1 - rdt_detectable(variables[[i]], parameters[[i]], timestep)
      }
    )
    p_mix <- mixing[[match_timestep(mixing_tt, timestep)]]
    p_captured_t <- p_captured[[match_timestep(p_captured_tt, timestep)]]
    n_species <- length(parameters[[1]]$species)
    species_eir <- t(vapply(
      seq_along(lagged_eir),
      function(i) {
        vnapply(
          lagged_eir[[i]],
          function(l) l$get(timestep - parameters[[i]]$de)
        )
      },
      numeric(n_species)
    ))
    if (n_species == 1) {
      species_eir <- t(species_eir)
    }
    inf <- vnapply(
      seq_along(lagged_infectivity),
      function(i) lagged_infectivity[[i]]$get(timestep - parameters[[i]]$delay_gam)
    )

    treated_coeff <- p_captured_t * rdt_negative * p_success

    eir <- vapply(
      seq(n_species),
      function(i) {
        mixed_eir <- t(species_eir[,i] * p_mix)
        eir_untreated <- mixed_eir * (1 - p_captured_t)
        eir_treated <- mixed_eir * treated_coeff
        rowSums(eir_untreated + eir_treated)
      },
      numeric(n_pops)
    )

    mixed_inf <- inf * p_mix
    inf_untreated <- mixed_inf * (1 - p_captured_t)
    inf_treated <- mixed_inf * treated_coeff

    inf <- rowSums(inf_untreated + inf_treated)

    list(eir = eir, inf = inf)
  }
}

# Estimates RDT prevalence from miscroscopy prevalence,
# Wu et al 2015 Nature paper
rdt_detectable <- function(variables, parameters, timestep) {
  microscopy_prev <- calculate_detected(
    variables$state,
    variables$birth,
    variables$id,
    parameters,
    timestep
  )$size() / parameters$human_population
  logit_prev <- log(microscopy_prev / (1 - microscopy_prev))
  logit_rdt <- parameters$rdt_intercept + parameters$rdt_coeff * logit_prev
  exp(logit_rdt) / (1 + exp(logit_rdt))
}

# return proportion of individuals detectable by microscopy
calculate_detected <- function(state, birth, immunity, parameters, timestep) {
  asymptomatic <- state$get_index_of('A')
  prob <- probability_of_detection(
    get_age(birth$get_values(asymptomatic), timestep),
    immunity$get_values(asymptomatic),
    parameters
  )
  asymptomatic_detected <- bitset_at(asymptomatic, bernoulli_multi_p(prob))
  clinically_detected <- state$get_index_of(c('Tr', 'D'))
  clinically_detected$or(asymptomatic_detected)
}
