#' @title Parameterise HRP2 antigen persistence
#' @description
#' Sets the weibull shape/scale parameters governing the hazard of losing
#' HRP2 positivity as a function of time since the triggering infection, the
#' probability that a new asymptomatic infection is HRP2 positive, and
#' whether treated infections clear HRP2 at a different rate to untreated/
#' asymptomatic infections. See \code{\link{get_parameters}} for details and
#' current (placeholder) defaults.
#' @param parameters model parameters
#' @param hrp2_shape weibull shape parameter for untreated/asymptomatic infections;
#' default = NULL (leave unchanged)
#' @param hrp2_scale weibull scale parameter (in timesteps) for untreated/asymptomatic
#' infections; default = NULL (leave unchanged)
#' @param hrp2_asymptomatic_prob probability that a new asymptomatic infection which is
#' detectable by light microscopy is also HRP2 positive; default = NULL (leave unchanged)
#' @param hrp2_treated_same_clearance if TRUE (the default), treated infections clear HRP2
#' at the same rate as untreated/asymptomatic infections (hrp2_shape/hrp2_scale); if FALSE,
#' treated infections instead use hrp2_shape_treated/hrp2_scale_treated;
#' default = NULL (leave unchanged)
#' @param hrp2_shape_treated weibull shape parameter for treated infections, only used
#' when hrp2_treated_same_clearance is FALSE; default = NULL (leave unchanged)
#' @param hrp2_scale_treated weibull scale parameter (in timesteps) for treated infections,
#' only used when hrp2_treated_same_clearance is FALSE; default = NULL (leave unchanged)
#' @export
set_hrp2_parameters <- function(
    parameters,
    hrp2_shape = NULL,
    hrp2_scale = NULL,
    hrp2_asymptomatic_prob = NULL,
    hrp2_treated_same_clearance = NULL,
    hrp2_shape_treated = NULL,
    hrp2_scale_treated = NULL
) {
  if (!is.null(hrp2_shape)) {
    parameters$hrp2_shape <- hrp2_shape
  }
  if (!is.null(hrp2_scale)) {
    parameters$hrp2_scale <- hrp2_scale
  }
  if (!is.null(hrp2_asymptomatic_prob)) {
    parameters$hrp2_asymptomatic_prob <- hrp2_asymptomatic_prob
  }
  if (!is.null(hrp2_treated_same_clearance)) {
    parameters$hrp2_treated_same_clearance <- hrp2_treated_same_clearance
  }
  if (!is.null(hrp2_shape_treated)) {
    parameters$hrp2_shape_treated <- hrp2_shape_treated
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
    if (parameters$parasite == "falciparum") {
      lm_detection_prob <- probability_of_detection(
        get_age(variables$birth$get_values(to_A), timestep),
        variables$id$get_values(to_A),
        parameters
      )
    } else {
      # p.v asymptomatic infections are always lm-detectable (see
      # create_prevalence_renderer)
      lm_detection_prob <- rep(1, to_A$size())
    }
    hrp2_prob <- lm_detection_prob * parameters$hrp2_asymptomatic_prob
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
#' individuals according to a weibull survival function of the time elapsed
#' since their triggering infection. The per-timestep clearance probability
#' is derived as the discrete hazard between consecutive timesteps:
#' 1 - S(x+1)/S(x).
#' @param variables a list of all of the model variables
#' @param timestep current timestep
#' @param positive bitset of hrp2-positive individuals in this group
#' @param shape weibull shape parameter for this group
#' @param scale weibull scale parameter for this group
#' @noRd
clear_hrp2_group <- function(variables, timestep, positive, shape, scale) {
  if (positive$size() == 0) {
    return()
  }
  x <- timestep - variables$hrp2_infection_time$get_values(positive)
  s_x <- weibull_survival(x, shape, scale)
  s_x1 <- weibull_survival(x + 1, shape, scale)
  # guard against 0/0 once survival has numerically underflowed to zero
  p_clear <- ifelse(s_x == 0, 1, 1 - s_x1 / s_x)
  cleared <- bitset_at(positive, bernoulli_multi_p(p_clear))
  if (cleared$size() > 0) {
    variables$hrp2_infection_time$queue_update(-1, cleared)
  }
}

#' @title HRP2 clearance process
#' @description
#' Each timestep, clears HRP2 positivity from individuals. If
#' parameters$hrp2_treated_same_clearance is TRUE, everyone clears at the
#' same rate (hrp2_shape/hrp2_scale). If FALSE, the treated group (see
#' update_hrp2_and_render_infections) clears at a separate rate
#' (hrp2_shape_treated/hrp2_scale_treated).
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
      clear_hrp2_group(variables, timestep, positive, parameters$hrp2_shape, parameters$hrp2_scale)
    } else {
      treated_positive <- variables$hrp2_treated$get_index_of('treated')$and(positive)
      untreated_positive <- variables$hrp2_treated$get_index_of('untreated')$and(positive)
      clear_hrp2_group(variables, timestep, untreated_positive, parameters$hrp2_shape, parameters$hrp2_scale)
      clear_hrp2_group(variables, timestep, treated_positive, parameters$hrp2_shape_treated, parameters$hrp2_scale_treated)
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
