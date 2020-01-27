
bind_process_to_default_model <- function(process, parameters) {
  states <- create_states()
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables)
  bind_process_to_model(process, individuals, states, variables)
}

mock_returns <- function(returns) {
  call <- 0
  function(...) {
    call <<- call + 1
    returns[[call]]
  }
}

expect_any <- function(X, FUN) {
  for (x in X) {
    if (FUN(x) == TRUE) {
      expect(TRUE, 'match found')
      return()
    }
  }
  expect(FALSE, 'No match')
}

mock_simulation_frame <- function(values) {
  list(
    get_state = function(individual, ...) {
      subset <- c()
      for (state in list(...)) {
        subset <- c(subset, values[[individual$name]][[state$name]])
      }
      subset
    },
    get_variable = function(individual, variable) {
      values[[individual$name]][[variable$name]]
    },
    get_constant = function(individual, variable) {
      values[[individual$name]][[variable$name]]
    }
  )
}
