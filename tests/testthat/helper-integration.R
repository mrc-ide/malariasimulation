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

mock_category <- function(...) {
  v <- individual::CategoricalVariable$new(...)
  list(
    get_index_of = v$get_index_of,
    queue_update = mockery::mock()
  )
}

mock_double <- function(...) {
  v <- individual::DoubleVariable$new(...)
  list(
    get_values = v$get_values,
    queue_update = mockery::mock()
  )
}

mock_render <- function(...) {
  v <- individual::Render$new(...)
  list(
    render = mockery::mock()
  )
}

mock_event <- function(event) {
  list(
    get_scheduled = function(...) event$get_scheduled(...),
    schedule = mockery::mock()
  )
}

expect_bitset_update <- function(mock, value, index, call = 1) {
  expect_equal(mockery::mock_args(mock)[[call]][[1]], value)
  expect_equal(mockery::mock_args(mock)[[call]][[2]]$to_vector(), index)
}

# Determine if range of vector is FP 0.
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
