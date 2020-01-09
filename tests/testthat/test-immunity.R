test_that('immunity_decay boosts after a delay', {
  level <- c(2.4, 1.2, 0)
  last_timestep <- c(11, 5, 1)
  timestep <- 15
  rate <- .5
  delay <- 4
  expect_equal(
    immunity_decay(level, last_timestep, timestep, rate, delay),
    c(3.4, .6, 0)
  )
})

test_that('immunity_decay works for untouched individuals', {
  level <- c(2.4, 1.2, 0)
  last_timestep <- c(11, 5, -1)
  timestep <- 15
  rate <- .5
  delay <- 4
  expect_equal(
    immunity_decay(level, last_timestep, timestep, rate, delay),
    c(3.4, .6, 0)
  )
})
