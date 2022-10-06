
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
