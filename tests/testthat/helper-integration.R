
expect_any <- function(X, FUN) {
  for (x in X) {
    if (FUN(x) == TRUE) {
      expect(TRUE, 'Match found')
      return()
    }
  }
  expect(FALSE, 'No match')
}

expect_none <- function(X, FUN) {
  for (x in X) {
    if (FUN(x) == TRUE) {
      expect(FALSE, 'Unexpected match found')
      return()
    }
  }
  expect(TRUE, 'No match found')
}

mock_api <- function(values, parameters = list(), timestep = 1) {
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
    queue_state_update = mockery::mock(),
    queue_variable_update = mockery::mock(),
    schedule = mockery::mock(),
    clear_schedule = mockery::mock(),
    get_scheduled = mockery::mock(),
    get_timestep = function() timestep,
    get_parameters = function() parameters,
    render = mockery::mock()
  )
}

expect_variable_update <- function(args, name, value, index) {
  expect_equal(args[[2]]$name, name)
  expect_equal(args[[3]], value)
  expect_equal(args[[4]], index)
}
