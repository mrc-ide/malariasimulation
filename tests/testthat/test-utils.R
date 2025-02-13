
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

test_that("merged_named_lists works with duplicate 'a's in lists", {
  x <- list(
    { tmp <- list(1, 2); names(tmp) <- c("a", "a"); tmp },
    { tmp <- list(3); names(tmp) <- "a"; tmp },
    { tmp <- list(3); names(tmp) <- "b"; tmp }
  )
  result <- merged_named_lists(x)
  expected <- list(a = 1, b = 3)
  expect_equal(result, expected)
})

test_that("merged_named_lists works with single list containing duplicates", {
  x <- list(
    { tmp <- list(1, 2, 3); names(tmp) <- c("a", "a", "b"); tmp }
  )
  result <- merged_named_lists(x)
  expected <- list(a = 1, b = 3)
  expect_equal(result, expected)
})

test_that("merged_named_lists works with mixed lists and top-level elements", {
  x <- list(
    list(a = 1, b = 2),
    list(a = 3),
    a = 4
  )
  result <- merged_named_lists(x)
  expected <- list(a = 1, b = 2)
  expect_equal(result, expected)
})

test_that("merged_named_lists works with multiple list arguments", {
  result <- merged_named_lists(list(a = 1, b = 2), list(a = 3), list(b = 3))
  expected <- list(a = 1, b = 2)
  expect_equal(result, expected)
})