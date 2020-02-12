test_that('carrying_capacity is calculated correctly', {
  parameters <- list(
    K0    = 10,
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    R_bar= 2,
    days_per_timestep = 1
  )
  expect_equal(
    carrying_capacity(100, parameters),
    5.63,
    tolerance=1e-2
  )
})

test_that('carrying_capacity is takes into account the timescale', {
  parameters <- list(
    K0    = 10,
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    R_bar= 2,
    days_per_timestep = 5
  )
  expect_equal(
    carrying_capacity(100, parameters),
    12.8,
    tolerance = 1e-1
  )
})

test_that('carrying_capacity cycles every year', {
  parameters <- list(
    K0    = 10,
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    R_bar= 2,
    days_per_timestep = 1
  )

  time_points <- c(1, 30, 160, 240, 365)
  for (t in time_points) {
    for (y in 1:3) {
      expect_equal(
        carrying_capacity(t, parameters),
        carrying_capacity(t + 365 * y, parameters),
        tolerance = 1e-1
      )
    }
  }
})
