
simulate_until_stable <- function(
  variables,
  events,
  processes,
  solvers,
  parameters,
  tolerance = 1e-1,
  years = 3,
  max_t = 500 * 365
  ) {
  t <- 1
  stop_fn <- stable_mean_EIR(
    variables,
    solvers,
    parameters,
    years=years,
    tolerance=tolerance
  )

  for (t in seq(max_t)) {
    for (process in processes) {
      individual:::execute_any_process(process, t)
    }
    for (event in events) {
      event$.process()
    }
    for (variable in variables) {
      variable$.update()
    }
    for (event in events) {
      event$.resize()
    }
    for (variable in variables) {
      variable$.resize()
    }
    for (event in events) {
      event$.tick()
    }
    if (stop_fn(t)) {
      return(t)
    }
  }
  stop(paste0('exiting at t=', max_t))
}

simulate_from_t <- function(
  processes,
  variables,
  events,
  start_t,
  timesteps
  ) {
  for (t in seq(start_t, start_t + timesteps)) {
    for (process in processes) {
      individual:::execute_any_process(process, t)
    }
    for (event in events) {
      event$.process()
    }
    for (variable in variables) {
      variable$.update()
    }
    for (event in events) {
      event$.resize()
    }
    for (variable in variables) {
      variable$.resize()
    }
    for (event in events) {
      event$.tick()
    }
  }
}


stable_mean_EIR <- function(
  variables,
  solvers,
  parameters,
  years=3,
  tolerance=1e-1
  ) {
  buffer = new.env()
  buffer$last_means <- rep(NA, years)
  buffer$last_year <- rep(NA, 365)
  function(t) {
    buffer$last_year[[(t %% 365) + 1]] <- overall_EIR(variables, solvers, parameters, t)
    if (t %% 365 == 0) {
      i <- ((floor(t / 365) - 1) %% years) + 1
      buffer$last_means[[i]] <- mean(buffer$last_year)
      if (any(is.na(buffer$last_means))) {
        return(FALSE)
      }
      order <- ((seq(years) + (i - 1)) %% years) + 1
      return(all(abs(diff(buffer$last_means[order])) < tolerance))
    }
    return(FALSE)
  }
}

overall_EIR <- function(variables, solvers, parameters, t) {
  sum(vnapply(
    seq_along(parameters$species),
    function (s) calculate_eir(s, solvers, variables, parameters, t)
  )) * 365 / parameters$human_population
}
