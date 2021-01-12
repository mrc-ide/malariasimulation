mock_double <- function(values) {
  list(
    get_values = mockery::mock(values),
    queue_update = mockery::mock()
  )
}

mock_category <- function(values) {
  list(
    get_index_of = function(categories) {
      b <- individual::Bitset$new(length(values))
      for (c in categories) {
        b$insert(which(values == c))
      }
      b
    },
    queue_update = mockery::mock()
  )
}

mock_renderer <- function() {
  list(
    render = mockery::mock()
  )
}

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

expect_variable_update <- function(args, name, value, index) {
  expect_equal(args[[2]]$name, name)
  expect_equal(args[[3]], value)
  expect_equal(args[[4]], index)
}

# Determine if range of vector is FP 0.
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
