test_that('competing_outcomes can handle no outcomes', {
  expect_equal(
    simulate_competing_outcomes(matrix(0, ncol=3, nrow=4)),
    rep(NA, 4)
  )
})

test_that('competing_outcomes can handle 1 outcomes', {
  mockery::stub(
    simulate_competing_outcomes,
    'runif',
    c(0, .1, 1, .2)
  )
  mockery::stub(
    simulate_competing_outcomes,
    'sample.int',
    mockery::mock(1)
  )
  expect_equal(
    simulate_competing_outcomes(
      matrix(
        c(
          0, .2, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0
        ),
        ncol=3,
        nrow=4
      )
    ),
    c(NA, 1, NA, NA)
  )
})

test_that('competing_outcomes can handle some outcomes', {
  mockery::stub(
    simulate_competing_outcomes,
    'runif',
    c(0, .1, 1, .2)
  )
  mockery::stub(
    simulate_competing_outcomes,
    'sample.int',
    mockery::mock(1, 3)
  )
  expect_equal(
    simulate_competing_outcomes(
      matrix(
        c(
          0, .2, 0, .4,
          0, 0, .3, 0,
          0, 0, 0, .1
        ),
        ncol=3,
        nrow=4
      )
    ),
    c(NA, 1, NA, 3)
  )
})
