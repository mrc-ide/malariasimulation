test_that('lag gives default when empty', {
  expect_equal(
    LaggedValue$new(12.5, 10)$get(1),
    10
  )
})

test_that('lag gives default before', {
  expect_equal(
    LaggedValue$new(12.5, 10)$get(1),
    10
  )
})

test_that('lag gives default before', {
  lagged <- LaggedValue$new(12.5, 10)
  lagged$save(2, 2)
  expect_equal(
    lagged$get(1),
    10
  )
})

test_that('lag gives default after', {
  lagged <- LaggedValue$new(12.5, 10)
  lagged$save(2, 2)
  expect_equal(
    lagged$get(3),
    10
  )
})

test_that('lag interpolates inbetween', {
  lagged <- LaggedValue$new(12.5, 10)
  lagged$save(2, 2)
  lagged$save(4, 4)
  expect_equal(
    lagged$get(3),
    3
  )
})

test_that('returns on the boundary', {
  lagged <- LaggedValue$new(12.5, 10)
  lagged$save(2, 2)
  expect_equal(
    lagged$get(2),
    2
  )
})

test_that('returns on the boundary (multiple points)', {
  lagged <- LaggedValue$new(12.5, 10)
  lagged$save(2, 2)
  lagged$save(4, 4)
  expect_equal(
    lagged$get(4),
    4
  )
})
