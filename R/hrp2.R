#' @title Parameterise HRP2 antigen persistence
#' @description
#' Sets the mean duration between an infection actually clearing and HRP2
#' becoming undetectable (a exponential/constant-hazard clearance process,
#' matching how infections are cleared elsewhere in the model - e.g.
#' dd/da/du/dt, see \code{\link{get_parameters}}), and the probability that
#' a new asymptomatic infection is HRP2 positive. See
#' \code{\link{get_parameters}} for details and current (placeholder)
#' defaults.
#' @param parameters model parameters
#' @param hrp2_scale mean duration between infection clearance and HRP2
#' clearance, in timesteps (the same rate is used for everyone - treated and
#' untreated individuals alike); default = NULL (leave unchanged)
#' @param hrp2_asymptomatic_prob probability that a new asymptomatic infection which is
#' detectable by light microscopy is also HRP2 positive; default = NULL (leave unchanged)
#' @export
set_hrp2_parameters <- function(
    parameters,
    hrp2_scale = NULL,
    hrp2_asymptomatic_prob = NULL
) {
  if (!is.null(hrp2_scale)) {
    parameters$hrp2_scale <- hrp2_scale
  }
  if (!is.null(hrp2_asymptomatic_prob)) {
    parameters$hrp2_asymptomatic_prob <- hrp2_asymptomatic_prob
  }
  parameters
}

#' @title Render new clinical infection counts and update HRP2 status
#' @description
#' Renders the number of new clinical infections this timestep, split into
#' treated and untreated, and marks new clinical infections as HRP2
#' positive. New asymptomatic (patent, non-clinical) infections (to_A) are
#' also marked HRP2 positive, but only with probability
#' `hrp2_asymptomatic_prob` * P(this asymptomatic infection is detectable by
#' light microscopy) - i.e. HRP2 is assumed to catch a subset of the
#' asymptomatic infections that light microscopy catches, not all of them.
#' Subpatent infections do not produce detectable HRP2 and are excluded.
#' Only to_D/treated are counted towards n_treated_infections/
#' n_untreated_infections: to_A never contributes to those counts, only to
#' HRP2 positivity.
#'
#' Marking someone HRP2 positive here does not, by itself, start any
#' clearance clock - see create_hrp2_clearance_process(), which only clears
#' HRP2 once the underlying infection has actually resolved (state == 'S'),
#' so reinfection while still HRP2 positive is a no-op.
#' @param variables a list of all of the model variables
#' @param renderer model render object
#' @param parameters model parameters
#' @param timestep current timestep
#' @param to_D bitset of humans newly moving to state D (untreated clinical)
#' @param treated bitset of humans newly moving to state Tr (treated clinical)
#' @param to_A bitset of humans newly moving to state A (asymptomatic, patent)
#' @noRd
update_hrp2_and_render_infections <- function(
    variables,
    renderer,
    parameters,
    timestep,
    to_D,
    treated,
    to_A
) {
  renderer$render('n_untreated_infections', to_D$size(), timestep)
  renderer$render('n_treated_infections', treated$size(), timestep)

  incidence_renderer(
    variables$birth,
    renderer,
    to_D,
    'untreated_infections_',
    parameters$hrp2_rendering_min_ages,
    parameters$hrp2_rendering_max_ages,
    timestep
  )
  incidence_renderer(
    variables$birth,
    renderer,
    treated,
    'treated_infections_',
    parameters$hrp2_rendering_min_ages,
    parameters$hrp2_rendering_max_ages,
    timestep
  )

  hrp2_positive_A <- to_A
  if (to_A$size() > 0) {
    hrp2_prob <- hrp2_asymptomatic_probability(
      get_age(variables$birth$get_values(to_A), timestep),
      if (parameters$parasite == "falciparum") variables$id$get_values(to_A) else NULL,
      parameters
    )
    hrp2_positive_A <- bitset_at(to_A, bernoulli_multi_p(hrp2_prob))
  }

  new_hrp2_positive <- to_D$copy()$or(treated)$or(hrp2_positive_A)
  if (new_hrp2_positive$size() > 0) {
    variables$hrp2$queue_update(1, new_hrp2_positive)
  }
}

#' @title HRP2 clearance process
#' @description
#' Each timestep, clears HRP2 positivity from individuals whose underlying
#' infection has actually resolved. HRP2-positive individuals still in an
#' infected state (D/A/U/Tr) are left untouched - they are still
#' parasitaemic, so HRP2 stays positive with no decay. Once an individual
#' reaches state 'S' (infection cleared), a constant per-timestep hazard of
#' 1/hrp2_scale is applied - the same exponential/memoryless convention used
#' for clearing infections elsewhere in the model (e.g.
#' progression_rates = 1/dd, 1/da, 1/du, 1/dt), rather than a hazard that
#' varies with time since clearance. The same hazard is used regardless of
#' whether the resolved infection was treated or not - treated individuals
#' clear HRP2 sooner only because their underlying infection (Tr -> S)
#' resolves faster than an untreated one (D -> A -> U -> S).
#' @param variables a list of all of the model variables
#' @param parameters model parameters
#' @noRd
create_hrp2_clearance_process <- function(variables, parameters) {
  function(timestep) {
    positive <- variables$hrp2$get_index_of(1)
    if (positive$size() == 0) {
      return()
    }
    cleared_infection <- variables$state$get_index_of('S')$and(positive)
    if (cleared_infection$size() == 0) {
      return()
    }
    p_clear <- rate_to_prob(1 / parameters$hrp2_scale)
    cleared <- bitset_at(
      cleared_infection,
      bernoulli_multi_p(rep(p_clear, cleared_infection$size()))
    )
    if (cleared$size() > 0) {
      variables$hrp2$queue_update(0, cleared)
    }
  }
}

#' @title Render the total number of HRP2 positive humans
#' @param renderer model render object
#' @param variables a list of all of the model variables
#' @noRd
create_hrp2_renderer_process <- function(renderer, variables) {
  function(timestep) {
    renderer$render(
      'n_hrp2_positive',
      variables$hrp2$get_index_of(1)$size(),
      timestep
    )
  }
}

#' @title Render the number of HRP2 positive humans by age band
#' @param variables a list of all of the model variables
#' @param parameters model parameters
#' @param renderer model render object
#' @noRd
create_hrp2_age_renderer_process <- function(variables, parameters, renderer) {
  function(timestep) {
    positive <- variables$hrp2$get_index_of(1)
    for (i in seq_along(parameters$hrp2_rendering_min_ages)) {
      lower <- parameters$hrp2_rendering_min_ages[[i]]
      upper <- parameters$hrp2_rendering_max_ages[[i]]
      in_age <- in_age_range(variables$birth, timestep, lower, upper)
      renderer$render(
        paste0('n_hrp2_positive_', lower, '_', upper),
        positive$copy()$and(in_age)$size(),
        timestep
      )
    }
  }
}

#' @title Probability that an asymptomatic infection is HRP2 positive
#' @description
#' P(HRP2 positive) = hrp2_asymptomatic_prob * P(this infection is
#' detectable by light microscopy). For p.f, LM detectability depends on age
#' and detection immunity (id); for p.v, asymptomatic infections are always
#' assumed LM-detectable (see create_prevalence_renderer()).
#' @param age age of each individual, in timesteps
#' @param id detection immunity of each individual (falciparum only; NULL for vivax)
#' @param parameters model parameters
#' @noRd
hrp2_asymptomatic_probability <- function(age, id, parameters) {
  if (parameters$parasite == "falciparum") {
    lm_detection_prob <- probability_of_detection(age, id, parameters)
  } else {
    lm_detection_prob <- rep(1, length(age))
  }
  lm_detection_prob * parameters$hrp2_asymptomatic_prob
}
