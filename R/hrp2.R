#' @title Parameterise HRP2 antigen persistence
#' @description
#' Sets the weibull shape/scale parameters governing the hazard of losing
#' HRP2 positivity as a function of time since the triggering clinical
#' infection. See \code{\link{get_parameters}} for details and current
#' (placeholder) defaults.
#' @param parameters model parameters
#' @param hrp2_shape weibull shape parameter; default = NULL (leave unchanged)
#' @param hrp2_scale weibull scale parameter, in timesteps; default = NULL (leave unchanged)
#' @export
set_hrp2_parameters <- function(
    parameters,
    hrp2_shape = NULL,
    hrp2_scale = NULL
) {
  if (!is.null(hrp2_shape)) {
    parameters$hrp2_shape <- hrp2_shape
  }
  if (!is.null(hrp2_scale)) {
    parameters$hrp2_scale <- hrp2_scale
  }
  parameters
}

#' @title Render new clinical infection counts and update HRP2 status
#' @description
#' Renders the number of new clinical infections this timestep, split into
#' treated and untreated, and (re)starts the HRP2 positivity clock for all
#' new clinical AND new asymptomatic (patent, non-clinical) infections - any
#' new blood-stage infection that isn't purely subpatent produces HRP2
#' antigen. Subpatent infections are not detectable and are excluded.
#' Only to_D/treated are counted towards n_treated_infections/
#' n_untreated_infections: to_A is used solely to (re)start the HRP2 clock.
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

  new_hrp2_positive <- to_D$copy()$or(treated)$or(to_A)
  if (new_hrp2_positive$size() > 0) {
    variables$hrp2_infection_time$queue_update(timestep, new_hrp2_positive)
  }
}

#' @title HRP2 clearance process
#' @description
#' Each timestep, clears HRP2 positivity from individuals according to a
#' weibull survival function of the time elapsed since their triggering
#' clinical infection. The per-timestep clearance probability is derived as
#' the discrete hazard between consecutive timesteps: 1 - S(x+1)/S(x).
#' @param variables a list of all of the model variables
#' @param parameters model parameters
#' @noRd
create_hrp2_clearance_process <- function(variables, parameters) {
  function(timestep) {
    positive <- variables$hrp2_infection_time$get_index_of(-1)$not(TRUE)
    if (positive$size() > 0) {
      x <- timestep - variables$hrp2_infection_time$get_values(positive)
      s_x <- weibull_survival(x, parameters$hrp2_shape, parameters$hrp2_scale)
      s_x1 <- weibull_survival(x + 1, parameters$hrp2_shape, parameters$hrp2_scale)
      # guard against 0/0 once survival has numerically underflowed to zero
      p_clear <- ifelse(s_x == 0, 1, 1 - s_x1 / s_x)
      cleared <- bitset_at(positive, bernoulli_multi_p(p_clear))
      if (cleared$size() > 0) {
        variables$hrp2_infection_time$queue_update(-1, cleared)
      }
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
