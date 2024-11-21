
test_that('time_cached calls the function the first time round', {
  fn <- mockery::mock(42)
  cached_fn <- time_cached(fn)
  expect_equal(cached_fn(timestep = 1), 42)
  mockery::expect_called(fn, 1)
})

test_that('time_cached returns the cached result the second time round', {
  fn <- mockery::mock(42)
  cached_fn <- time_cached(fn)
  expect_equal(cached_fn(timestep = 1), 42)
  expect_equal(cached_fn(timestep = 1), 42)
  mockery::expect_called(fn, 1)
})

test_that('time_cached calls again for a different timestep', {
  fn <- mockery::mock(42, 43)
  cached_fn <- time_cached(fn)
  expect_equal(cached_fn(timestep = 1), 42)
  expect_equal(cached_fn(timestep = 2), 43)
  mockery::expect_called(fn, 2)
})

test_that("bitset_index works", {
  a <- individual::Bitset$new(10)$insert(c(3,5,7,9))
  b <- individual::Bitset$new(10)$insert(c(2,4,5,8,9))
  expect_equal(bitset_index(a, b), c(2,4))

  a <- individual::Bitset$new(10)
  b <- individual::Bitset$new(10)$insert(c(2,4,5,8,9))
  expect_equal(length(bitset_index(a, b)), 0)

  a <- individual::Bitset$new(10)$insert(c(3,5,7,9))
  b <- individual::Bitset$new(10)
  expect_equal(length(bitset_index(a, b)), 0)
})

test_that("bitset_index errors if size does not match", {
  a <- individual::Bitset$new(10)$insert(c(3,5,7,9))
  b <- individual::Bitset$new(20)$insert(c(2,4,5,8,9))
  expect_error(bitset_index(a, b), "Incompatible bitmap sizes")
})
