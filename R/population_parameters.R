#' @title Parameterise variable deathrates
#'
#' @param parameters the model parameters
#' @param agegroups vector of agegroups (in timesteps)
#' @param timesteps vector of timesteps for each change in demography
#' @param deathrates matrix of deathrates per age group per timestep.
#' Rows are timesteps from the `timesteps` param. Columns are the age groups
#' from the `agegroups` param.
#' @export
set_demography <- function(
  parameters,
  agegroups,
  timesteps,
  deathrates
  ) {
  
  stopifnot(all(timesteps >= 0))
  if(min(timesteps) != 0){
    stop("Must include the baseline demography (timesteps must include 0),
         when setting a custom demography")
  }
  stopifnot(all(agegroups > 0))
  stopifnot(all(deathrates > 0 & deathrates < 1))
  stopifnot(length(agegroups) == ncol(deathrates))
  stopifnot(length(timesteps) == nrow(deathrates))
  stopifnot(!is.unsorted(timesteps, strictly = TRUE))

  parameters$custom_demography <- TRUE
  parameters$deathrate_agegroups <- agegroups
  parameters$deathrate_timesteps <- timesteps
  parameters$deathrates <- deathrates

  parameters
}

get_equilibrium_population <- function(age_high, birthrate, deathrates) {
  n_age <- length(age_high)
  # r[i] can be thought of as the rate of ageing in this age group, i.e.
  # 1/r[i] is the duration of this group
  age_width <- diff(c(0, age_high))
  aging_rate <- 1 / age_width
  aging_rate[[n_age]] <- 0

  pop <- rep(0, n_age)
  for (i in seq(n_age)) {
    # calculate poportions in each age group
    # birth rate inflow / (aging out + deathrate)
    if (i == 1) {
      pop[i] <- birthrate / (aging_rate[[i]] + deathrates[[i]])
    } else {
      pop[i] <- pop[[i-1]] * aging_rate[[i-1]] / (aging_rate[i] + deathrates[[i]])
    }
  }
  pop
}

#' @title Calculate the birthrate for a population in equilibrium
#'
#' @param pops a vector of populations
#' @param age_high a vector of age groups
#' @param deathrates vector of deathrates for each age group
#' @importFrom stats uniroot
find_birthrates <- function(pops, age_high, deathrates) {
  vnapply(
    pops,
    function(pop) {
      birth_to_pop <- function(b) {
        pop - sum(get_equilibrium_population(age_high, b, deathrates))
      }
      uniroot(birth_to_pop, c(0, pop))$root
    }
  )
}

get_human_population <- function(parameters, timestep) {
  last_pop <- match_last_timestep(parameters$human_population_timesteps, timestep)
  parameters$human_population[last_pop]
}
