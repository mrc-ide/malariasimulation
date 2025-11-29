#' @title Parameterise non-malarial fevers
#' @description
#' Set age specific daily rates of non-malarial fever.
#'
#' @param parameters model parameters
#' @param ages vector of upper age bounds (in timesteps) for each fever rate
#' @param rates vector of rates of non-malarial fever per day
#' @return updated parameter list
#' @export
set_nmf <- function(parameters, ages, rates){
  stopifnot(length(ages) == length(rates))
  stopifnot(all(rates >= 0))
  stopifnot(all(diff(c(0, ages)) > 0))
  parameters$nmf_ages <- ages
  parameters$nmf_rates <- rates
  parameters
}

#' @title Non malarial fever process
#' @description
#' Samples individuals experiencing non-malarial fever.
#'
#' @param variables list of model variables
#' @param parameters model parameters
#' @param renderer model renderer
#' @param nmf bitset used to store individuals with fever
#' @return a function called at each timestep
#' @noRd
create_nmf_process <- function(variables, parameters, renderer, nmf){
  renderer$set_default('n_nmf', 0)
  function(timestep){
    nmf$clear()
    pop <- get_human_population(parameters, timestep)
    if(length(parameters$nmf_rates) == 0 || pop == 0){
      return()
    }
    age <- get_age(variables$birth$get_values(), timestep)
    groups <- .bincode(age, c(0, parameters$nmf_ages))
    rates <- rep(0, pop)
    in_range <- !is.na(groups)
    rates[in_range] <- parameters$nmf_rates[groups[in_range]]
    nmf$insert(bernoulli_multi_p(rates))
    renderer$render('n_nmf', nmf$size(), timestep)
  }
}
