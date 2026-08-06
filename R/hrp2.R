#' @title Parameterise HRP2 antigen persistence
#' @description
#' Sets the mean duration of HRP2 positivity (an exponential/constant-hazard
#' clearance process, matching how infections are cleared elsewhere in the
#' model - e.g. dd/da/du/dt, see \code{\link{get_parameters}}), the
#' probability that a new asymptomatic infection is HRP2 positive, and
#' whether treated infections clear HRP2 at a different rate to untreated/
#' asymptomatic infections. See \code{\link{get_parameters}} for details and
#' current (placeholder) defaults.
#' @param parameters model parameters
#' @param hrp2_scale mean duration of HRP2 positivity for untreated/asymptomatic
#' infections, in timesteps; default = NULL (leave unchanged)
#' @param hrp2_asymptomatic_prob probability that a new asymptomatic infection which is
#' detectable by light microscopy is also HRP2 positive; default = NULL (leave unchanged)
#' @param hrp2_treated_same_clearance if TRUE (the default), treated infections clear HRP2
#' at the same rate as untreated/asymptomatic infections (hrp2_scale); if FALSE,
#' treated infections instead use hrp2_scale_treated;
#' default = NULL (leave unchanged)
#' @param hrp2_scale_treated mean duration of HRP2 positivity for treated infections, in
#' timesteps, only used when hrp2_treated_same_clearance is FALSE;
#' default = NULL (leave unchanged)
#' @export
set_hrp2_parameters <- function(
    parameters,
    hrp2_scale = NULL,
    hrp2_asymptomatic_prob = NULL,
    hrp2_treated_same_clearance = NULL,
    hrp2_scale_treated = NULL
) {
  if (!is.null(hrp2_scale)) {
    parameters$hrp2_scale <- hrp2_scale
  }
  if (!is.null(hrp2_asymptomatic_prob)) {
    parameters$hrp2_asymptomatic_prob <- hrp2_asymptomatic_prob
  }
  if (!is.null(hrp2_treated_same_clearance)) {
    parameters$hrp2_treated_same_clearance <- hrp2_treated_same_clearance
  }
  if (!is.null(hrp2_scale_treated)) {
    parameters$hrp2_scale_treated <- hrp2_scale_treated
  }
  parameters
}

#' @title Render new clinical infection counts and update HRP2 status
#' @description
#' Renders the number of new clinical infections this timestep, split into
#' treated and untreated, and (re)starts the HRP2 positivity clock for all
#' new clinical infections. New asymptomatic (patent, non-clinical)
#' infections (to_A) also (re)start the clock, but only with probability
#' `hrp2_asymptomatic_prob` * P(this asymptomatic infection is detectable by
#' light microscopy) - i.e. HRP2 is assumed to catch a subset of the
#' asymptomatic infections that light microscopy catches, not all of them.
#' Subpatent infections do not produce detectable HRP2 and are excluded.
#' Only to_D/treated are counted towards n_treated_infections/
#' n_untreated_infections: to_A never contributes to those counts, only to
#' the HRP2 clock.
#'
#' Whichever of the treated/untreated groups triggers a given individual's
#' HRP2 positivity is recorded (hrp2_treated), so create_hrp2_clearance_process
#' can later apply a distinct clearance hazard per group if
#' parameters$hrp2_treated_same_clearance is FALSE. Asymptomatic infections
#' (to_A) are grouped with untreated (to_D) for this purpose.
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

  untreated_hrp2_positive <- to_D$copy()$or(hrp2_positive_A)
  if (untreated_hrp2_positive$size() > 0) {
    variables$hrp2_infection_time$queue_update(timestep, untreated_hrp2_positive)
    variables$hrp2_treated$queue_update('untreated', untreated_hrp2_positive)
  }
  if (treated$size() > 0) {
    variables$hrp2_infection_time$queue_update(timestep, treated)
    variables$hrp2_treated$queue_update('treated', treated)
  }
}

#' @title Clear HRP2 positivity from a group of individuals
#' @description
#' Clears HRP2 positivity from the given (already hrp2-positive) subset of
#' individuals via a constant per-timestep hazard of 1/scale - the same
#' exponential/memoryless convention used for clearing infections elsewhere
#' in the model (e.g. progression_rates = 1/dd, 1/da, 1/du, 1/dt), rather
#' than a hazard that varies with time since infection.
#' @param variables a list of all of the model variables
#' @param positive bitset of hrp2-positive individuals in this group
#' @param scale mean duration of HRP2 positivity for this group, in timesteps
#' @noRd
clear_hrp2_group <- function(variables, positive, scale) {
  if (positive$size() == 0) {
    return()
  }
  p_clear <- rate_to_prob(1 / scale)
  cleared <- bitset_at(positive, bernoulli_multi_p(rep(p_clear, positive$size())))
  if (cleared$size() > 0) {
    variables$hrp2_infection_time$queue_update(-1, cleared)
  }
}

#' @title HRP2 clearance process
#' @description
#' Each timestep, clears HRP2 positivity from individuals. If
#' parameters$hrp2_treated_same_clearance is TRUE, everyone clears at the
#' same rate (hrp2_scale). If FALSE, the treated group (see
#' update_hrp2_and_render_infections) clears at a separate rate
#' (hrp2_scale_treated).
#' @param variables a list of all of the model variables
#' @param parameters model parameters
#' @noRd
create_hrp2_clearance_process <- function(variables, parameters) {
  function(timestep) {
    positive <- variables$hrp2_infection_time$get_index_of(-1)$not(TRUE)
    if (positive$size() == 0) {
      return()
    }
    if (parameters$hrp2_treated_same_clearance) {
      clear_hrp2_group(variables, positive, parameters$hrp2_scale)
    } else {
      treated_positive <- variables$hrp2_treated$get_index_of('treated')$and(positive)
      untreated_positive <- variables$hrp2_treated$get_index_of('untreated')$and(positive)
      clear_hrp2_group(variables, untreated_positive, parameters$hrp2_scale)
      clear_hrp2_group(variables, treated_positive, parameters$hrp2_scale_treated)
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
      variables$hrp2_infection_time$get_index_of(-1)$not(TRUE)$size(),
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
    positive <- variables$hrp2_infection_time$get_index_of(-1)$not(TRUE)
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
