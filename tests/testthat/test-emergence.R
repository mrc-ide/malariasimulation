test_that('carrying_capacity is calculated correctly', {
  parameters <- list(
    model_seasonality = TRUE,
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    days_per_timestep = 1
  )
  expect_equal(
    carrying_capacity(100, parameters, 10, 2),
    5.63,
    tolerance=1e-2
  )
})

test_that('carrying_capacity is takes into account the timescale', {
  parameters <- list(
    model_seasonality = TRUE,
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    days_per_timestep = 5
  )
  expect_equal(
    carrying_capacity(100, parameters, 10, 2),
    12.8,
    tolerance = 1e-1
  )
})

test_that('carrying_capacity cycles every year', {
  parameters <- list(
    model_seasonality = TRUE,
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    days_per_timestep = 1
  )

  time_points <- c(1, 30, 160, 240, 365)
  for (t in time_points) {
    for (y in 1:3) {
      expect_equal(
        carrying_capacity(t, parameters, 10, 2),
        carrying_capacity(t + 365 * y, parameters, 10, 2),
        tolerance = 1e-1
      )
    }
  }
})

test_that('carrying_capacity can avoid seasonality', {
  parameters <- list(
    model_seasonality = FALSE,
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    days_per_timestep = 1
  )

  time_points <- seq_len(5)
  for (t in time_points) {
    for (y in 1:3) {
      expect_equal(
        100,
        carrying_capacity(t + 365 * y, parameters, 100, 2),
        tolerance = 1e-1
      )
    }
  }
})
