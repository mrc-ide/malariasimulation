
simulate_until_stable <- function(
  variables,
  events,
  processes,
  solvers,
  parameters,
  tolerance = 1e-2,
  max_t = 500 * 365
  ) {
  t <- 1

  eir <- sum(vnapply(
    seq_along(parameters$species),
    function (s) calculate_eir(s, solvers, variables, parameters, t)
  )) / parameters$human_population
  eir_t <- NULL

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
    eir_t <- c(
      eir_t,
      sum(vnapply(
        seq_along(parameters$species),
        function (s) calculate_eir(s, solvers, variables, parameters, t)
      ))
    )
    if (t %% 365 == 0) {
      new_eir <- mean(eir_t) / parameters$human_population
      print(abs(new_eir - eir))
      if (abs(new_eir - eir) < tolerance) {
        return()
      }
      eir_t <- NULL
      eir <- new_eir
    }
  }
  stop(paste0('exiting at t=', max_t))
}
