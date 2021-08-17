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

mock_method <- function(class, method, mock) {
  MockClass <- R6::R6Class(
    paste0('Mock', class$classname),
    inherit = class
  )
  MockClass$set('public', method, function(...) mock(...))
  MockClass$set('public', paste0(method, '_mock'), function() mock)
  MockClass
}

mock_category <- function(...) {
  class <- mock_method(
    individual::CategoricalVariable,
    'queue_update',
    mockery::mock()
  )
  class$new(...)
}

mock_double <- function(...) {
  class <- mock_method(
    individual::DoubleVariable,
    'queue_update',
    mockery::mock()
  )
  class$new(...)
}

mock_render <- function(...) {
  class <- mock_method(
    individual::Render,
    'render',
    mockery::mock()
  )
  class$new(...)
}

mock_integer <- function(...) {
  class <- mock_method(
    individual::IntegerVariable,
    'queue_update',
    mockery::mock()
  )
  class$new(...)
}

mock_event <- function(event) {
  list(
    get_scheduled = function(...) event$get_scheduled(...),
    schedule = mockery::mock(),
    clear_schedule = mockery::mock()
  )
}

expect_bitset_update <- function(mock, value, index, call = 1) {
  expect_equal(mockery::mock_args(mock)[[call]][[1]], value)
  expect_equal(mockery::mock_args(mock)[[call]][[2]]$to_vector(), index)
}

expect_bitset_schedule <- function(mock, target, delay, call = 1) {
  expect_equal(mockery::mock_args(mock)[[call]][[1]]$to_vector(), target)
  expect_equal(mockery::mock_args(mock)[[call]][[2]], delay)
}

# Determine if range of vector is FP 0.
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
