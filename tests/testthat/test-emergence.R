test_that('carrying_capacity is calculated correctly', {
  expect_equal(
    carrying_capacity(
      100,
      TRUE,
      2,
      c(.3, .6, .9),
      c(.1, .4, .7),
      10,
      2,
      .001
    ),
    5.63,
    tolerance=1e-2
  )
})

test_that('carrying_capacity cycles every year', {
  time_points <- c(1, 30, 160, 240, 365)
  for (t in time_points) {
    for (y in 1:3) {
      expect_equal(
        carrying_capacity(
          t,
          TRUE,
          2,
          c(.3, .6, .9),
          c(.1, .4, .7),
          10,
          2,
          .001
        ),
        carrying_capacity(
          t + 365 * y,
          TRUE,
          2,
          c(.3, .6, .9),
          c(.1, .4, .7),
          10,
          2,
          .001
        ),
        tolerance = 1e-1
      )
    }
  }
})

test_that('carrying_capacity can avoid seasonality', {
  time_points <- seq_len(5)
  for (t in time_points) {
    for (y in 1:3) {
      expect_equal(
        100,
        carrying_capacity(
          t + 365 * y,
          FALSE,
          2,
          c(.3, .6, .9),
          c(.1, .4, .7),
          100,
          2,
          .001
        ),
        tolerance = 1e-1
      )
    }
  }
})
