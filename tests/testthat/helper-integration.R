
bind_process_to_default_model <- function(process, parameters) {
  states <- create_states()
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables)
  bind_process_to_model(process, individuals, states, variables)
}

mock_returns <- function(returns) {
  call <- 0
  function() {
    call <<- call + 1
    returns[[call]]
  }
}

expect_has_update <- function(update_list, update) {
  if(is.null(update_list)) {
    stop('Update list is NULL')
  }

  if(!is.list(update_list)) {
    update_list <- list(update_list)
  }

  has <- FALSE

  for (u in update_list) {
    if (updates_equal(update, u)) {
      has <- TRUE
    }
  }

  if (!has) {
    stop('update does not exist')
  }
}

updates_equal <- function(self, other) {
  if (inherits(self, 'StateUpdate')) {
    return(all(
      inherits(other, 'StateUpdate'),
      self$individual$name == other$individual$name,
      self$state$name == other$state$name,
      all.equal(self$index, other$index) == TRUE
    ))
  }
  all(
    inherits(other, 'VariableUpdate'),
    self$individual$name == other$individual$name,
    self$variable$name == other$variable$name,
    all.equal(self$value, other$value) == TRUE,
    all.equal(self$index, other$index) == TRUE
  )
}

mock_simulation_frame <- function(values) {
  list(
    get_state = function(individual, ...) {
      unlist(lapply(
        list(...),
        function(state) values[[individual$name]][[state$name]]
      ))
    },
    get_variable = function(individual, variable) {
      values[[individual$name]][[variable$name]]
    }
  )
}
