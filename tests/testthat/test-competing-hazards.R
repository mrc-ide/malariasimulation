# ==============================
# Test the CompetingHazard class
# ==============================

test_that('hazard resolves two normalised outcomes when all events occur', {

    size <- 4
    outcome_1_process <- mockery::mock()
    outcome_1 <- CompetingOutcome$new(
      targeted_process = outcome_1_process,
      size = size
    )
    outcome_2_process <- mockery::mock()
    outcome_2 <- CompetingOutcome$new(
      targeted_process = outcome_2_process,
      size = size
    )

    hazard <- CompetingHazard$new(
      outcomes = list(outcome_1, outcome_2),
      rng = mockery::mock(
        c(0, 0, 0, 0), # all events occur
        c(.05, .3, .2, .5) # event_rng
      )
    )

    outcome_1$set_rates(c(0.1, 0.2, 0.3, 0.4))
    outcome_2$set_rates(c(0.9, 0.8, 0.7, 0.6))

    hazard$resolve(0)

    mockery::expect_args(
      outcome_1_process,
      1,
      0,
      individual::Bitset$new(size)$insert(c(1, 3))
    )
    mockery::expect_args(
      outcome_2_process,
      1,
      0,
      individual::Bitset$new(size)$insert(c(2, 4))
    )
})

test_that('hazard resolves two unnormalised outcomes when all events occur', {

  size <- 4
  outcome_1_process <- mockery::mock()
  outcome_1 <- CompetingOutcome$new(
    targeted_process = outcome_1_process,
    size = size
  )

  outcome_2_process <- mockery::mock()
  outcome_2 <- CompetingOutcome$new(
    targeted_process = outcome_2_process,
    size = size
  )

  hazard <- CompetingHazard$new(
    outcomes = list(outcome_1, outcome_2),
    rng = mockery::mock(
      c(0, 0, 0, 0), # all events occur
      c(.05, .3, .2, .5) # event_rng
    )
  )

  outcome_1$set_rates(c(0.1, 0.2, 0.3, 0.4) / 2)
  outcome_2$set_rates(c(0.9, 0.8, 0.7, 0.6) / 2)

  hazard$resolve(0)

  mockery::expect_args(
    outcome_1_process,
    1,
    0,
    individual::Bitset$new(size)$insert(c(1, 3))
  )
  mockery::expect_args(
    outcome_2_process,
    1,
    0,
    individual::Bitset$new(size)$insert(c(2, 4))
  )
})

test_that('hazard resolves two outcomes when some events occur', {

  size <- 4
  outcome_1_process <- mockery::mock()
  outcome_1 <- CompetingOutcome$new(
    targeted_process = outcome_1_process,
    size = size
  )

  outcome_2_process <- mockery::mock()
  outcome_2 <- CompetingOutcome$new(
    targeted_process = outcome_2_process,
    size = size
  )

  hazard <- CompetingHazard$new(
    outcomes = list(outcome_1, outcome_2),
    rng = mockery::mock(
      c(0, 0, 1, 1), # some events occur
      c(.05, .3, .2, .5), # event_rng
    )
  )

  outcome_1$set_rates(c(0.1, 0.2, 0.3, 0.4) / 2)
  outcome_2$set_rates(c(0.9, 0.8, 0.7, 0.6) / 2)

  hazard$resolve(0)

  mockery::expect_args(
    outcome_1_process,
    1,
    0,
    individual::Bitset$new(size)$insert(1)
  )
  mockery::expect_args(
    outcome_2_process,
    1,
    0,
    individual::Bitset$new(size)$insert(2)
  )
})


test_that('hazard resolves when no rates are set', {

  size <- 4
  outcome_1_process <- mockery::mock()
  outcome_1 <- CompetingOutcome$new(
    targeted_process = outcome_1_process,
    size = size
  )

  outcome_2_process <- mockery::mock()
  outcome_2 <- CompetingOutcome$new(
    targeted_process = outcome_2_process,
    size = size
  )

  hazard <- CompetingHazard$new(
    outcomes = list(outcome_1, outcome_2)
  )

  hazard$resolve(0)

  mockery::expect_called(
    outcome_1_process,
    0
  )
  mockery::expect_called(
    outcome_2_process,
    0
  )
})

test_that('hazard resolves when rates are set to one', {

  size <- 4
  outcome_1_process <- mockery::mock()
  outcome_1 <- CompetingOutcome$new(
    targeted_process = outcome_1_process,
    size = size
  )

  outcome_2_process <- mockery::mock()
  outcome_2 <- CompetingOutcome$new(
    targeted_process = outcome_2_process,
    size = size
  )

  hazard <- CompetingHazard$new(
    outcomes = list(outcome_1, outcome_2)
  )

  outcome_1$set_rates(rep(1, size))
  hazard$resolve(0)

  mockery::expect_args(
    outcome_1_process,
    1,
    0,
    individual::Bitset$new(size)$insert(c(1, 2, 3, 4))
  )
  mockery::expect_called(
    outcome_2_process,
    0
  )
})

test_that('hazard resolves three outcomes', {

  size <- 4
  outcome_1_process <- mockery::mock()
  outcome_1 <- CompetingOutcome$new(
    targeted_process = outcome_1_process,
    size = size
  )

  outcome_2_process <- mockery::mock()
  outcome_2 <- CompetingOutcome$new(
    targeted_process = outcome_2_process,
    size = size
  )

  outcome_3_process <- mockery::mock()
  outcome_3 <- CompetingOutcome$new(
    targeted_process = outcome_3_process,
    size = size
  )

  hazard <- CompetingHazard$new(
    outcomes = list(outcome_1, outcome_2, outcome_3),
    rng = mockery::mock(
      c(0, 0, 0, 0), # all events occur
      c(.04, .3, .14, .6), # event_rng
    )
  )

  outcome_1$set_rates(c(0.1, 0.2, 0.3, 0.4) / 2)
  outcome_2$set_rates(c(0.9, 0.8, 0.7, 0.6) / 2)
  outcome_3$set_rates(c(0.5, 0.5, 0.5, 0.5))

  hazard$resolve(0)

  mockery::expect_args(
    outcome_1_process,
    1,
    0,
    individual::Bitset$new(size)$insert(c(1, 3))
  )
  mockery::expect_args(
    outcome_2_process,
    1,
    0,
    individual::Bitset$new(size)$insert(2)
  )
  mockery::expect_args(
    outcome_3_process,
    1,
    0,
    individual::Bitset$new(size)$insert(4)
  )
})

