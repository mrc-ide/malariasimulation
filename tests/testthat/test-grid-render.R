test_that('grid_count produces all age groups between 0 and 100', {
  timestep <- 5 * 365
  birth <- individual::IntegerVariable$new(c(1, 2, 3, 4) * 365)
  selected <- individual::Bitset$new(4)$not()
  expect_equal(length(grid_count(birth, selected, timestep)), 101)
})

test_that('grid_count works without selection', {
  timestep <- 5 * 365
  birth <- individual::IntegerVariable$new(c(1, 2, 3, 4) * 365)
  selected <- NULL
  expected <- rep(0, 101)
  expected[seq(4) + 1] <- 1
  expect_equal(grid_count(birth, selected, timestep), expected)
})

test_that('grid_count counts internally correctly', {
  timestep <- 5 * 365
  birth <- individual::IntegerVariable$new(timestep - c(1, 2, 3, 4) * 365 - 1)
  selected <- individual::Bitset$new(4)$not()
  expected <- rep(0, 101)
  expected[seq(4) + 1] <- 1
  expect_equal(grid_count(birth, selected, timestep), expected)
})

test_that('grid_count counts at the boundaries correctly', {
  timestep <- 5 * 365
  birth <- individual::IntegerVariable$new(timestep - c(1, 2, 3, 4) * 365)
  selected <- individual::Bitset$new(4)$not()
  expected <- rep(0, 101)
  expected[seq(4) + 1] <- 1
  expect_equal(grid_count(birth, selected, timestep), expected)
})

test_that('grid_count can select subset', {
  timestep <- 5 * 365
  birth <- individual::IntegerVariable$new(timestep - c(1, 2, 3, 4) * 365 - 1)
  selected <- individual::Bitset$new(4)$insert(c(2, 4))
  expected <- rep(0, 101)
  expected[c(3, 5)] <- 1
  expect_equal(grid_count(birth, selected, timestep), expected)
})

test_that('grid_count ignores outside of the grid', {
  timestep <- 5 * 365
  birth <- individual::IntegerVariable$new(timestep - c(1, 101, 3, 150) * 365)
  selected <- individual::Bitset$new(4)$not()
  expected <- rep(0, 101)
  expected[c(2, 4)] <- 1
  expect_equal(grid_count(birth, selected, timestep), expected)
})
