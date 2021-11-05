CONTINGENCY <- 1.1

#' @title Parameterise variable deathrates
#'
#' @param parameters the model parameters
#' @param agegroups vector of agegroups (in timesteps)
#' @param timesteps vector of timesteps for each change in demography
#' @param birthrates vector of new birthrates (number of individuals born per
#' timestep) for each element of `timesteps`
#' @param deathrates matrix of deathrates per age group per timestep.
#' Rows are timesteps from the `timesteps` param. Columns are the age groups
#' from the `agegroups` param.
#' @export
set_demography <- function(
  parameters,
  agegroups,
  timesteps,
  birthrates,
  deathrates
  ) {
  stopifnot(all(agegroups > 0))
  stopifnot(all(timesteps > 0))
  stopifnot(all(birthrates > 0))
  stopifnot(length(birthrates) == length(timesteps))
  stopifnot(all(deathrates > 0 & deathrates < 1))
  stopifnot(length(agegroups) == ncol(deathrates))
  stopifnot(length(timesteps) == nrow(deathrates))
  parameters$custom_demography <- TRUE
  parameters$deathrate_agegroups <- agegroups
  parameters$deathrate_timesteps <- timesteps
  parameters$deathrates <- deathrates
  parameters$birthrate_timesteps <- timesteps
  parameters$birthrates <- birthrates
  n_age <- length(agegroups)

  # set a max population
  populations <- vapply(
    seq(timesteps),
    function(timestep) {
      get_equilibrium_population(
        agegroups,
        birthrates[[timestep]],
        parameters$deathrates[timestep,]
      )
    },
    numeric(n_age)
  )

  deathrates <- vnapply(
    seq(timesteps),
    function(timestep) {
      sum(
         parameters$deathrates[timestep,] * populations[timestep,]
      )
    }
  )
  populations <- colSums(populations)
  parameters$max_human_population <- ceiling(
    max(populations * (1 + deathrates)) * CONTINGENCY
  )

  parameters
}

#' @title Parameterise the population size
#'
#' @description will set the initial population and estimate a max population
#' from parameterised birth/death rates
#'
#' @param parameters the model parameters
#' @param population initial population
set_initial_population <- function(parameters, population) {
  parameters$init_human_population <- population
  parameters$custom_demography <- FALSE
  deathrate <- 1 / parameters$average_age
  parameters$max_human_population <- ceiling(
    population * (1 + deathrate) * CONTINGENCY
  )
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

get_birthrate <- function(parameters, timestep) {
  if (!parameters$custom_demography) {
    return(1 / parameters$average_age * parameters$init_human_population)
  }
  print(which(timestep >= parameters$birthrate_timesteps))
  last_birthrate <- which(timestep >= parameters$birthrate_timesteps)[[-1]]
  parameters$birthrates[last_birthrate]
}
