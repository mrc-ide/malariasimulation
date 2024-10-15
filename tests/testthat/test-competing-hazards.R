# ==============================
# Test the CompetingHazard class
# ==============================


test_that("hazard resolves two disjoint outcomes", {
  size <- 4
  population <- individual::Bitset$new(size)$not()

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
    size = size,
    outcomes = list(outcome_1, outcome_2),
    rng = mockery::mock(c(.05, .3, .2, .5))
  )
  
  outcome_1$set_rates(population, c(10, 0, 10, 0))
  outcome_2$set_rates(population, c(0, 10, 0, 10))

  hazard$resolve(0)

  mockery::expect_args(outcome_1_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(1, 3)))

  mockery::expect_args(outcome_2_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(2, 4)))
  
})

test_that("hazard resolves two competing outcomes", {
  size <- 4
  population <- individual::Bitset$new(size)$not()

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
    size = size,
    outcomes = list(outcome_1, outcome_2),
    rng = mockery::mock(c(.7, .3, .2, .6))
  )
  
  outcome_1$set_rates(population, c(5, 5, 5, 5))
  outcome_2$set_rates(population, c(5, 5, 5, 5))
  
  hazard$resolve(0)
   
  mockery::expect_args(outcome_1_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(2, 3)))
  mockery::expect_args(outcome_2_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(1, 4)))
})

test_that("hazard may resolve to neither outcome", {
  size <- 4
  population <- individual::Bitset$new(size)$not()

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
    size = size,
    outcomes = list(outcome_1, outcome_2),
    rng = mockery::mock(c(.8, .4, .2, .6))
  )
  
  outcome_1$set_rates(population, prob_to_rate(rep(0.5, size)))
  outcome_2$set_rates(population, prob_to_rate(rep(0.5, size)))
  
  hazard$resolve(0)
   
  mockery::expect_args(outcome_1_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(3)))
  mockery::expect_args(outcome_2_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(2, 4)))
})

test_that("outcomes can define a partial set of rates", {
  size <- 4
  population <- individual::Bitset$new(size)$not()

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
    size = size,
    outcomes = list(outcome_1, outcome_2),

    # Only individuals 1, 3 and 4 get sampled
    rng = mockery::mock(c(.2, .3, .6))
  )
  
  outcome_1$set_rates(individual::Bitset$new(size)$insert(c(1,3)), c(5, 5))
  outcome_2$set_rates(individual::Bitset$new(size)$insert(c(1,4)), c(5, 5))
  
  hazard$resolve(0)
   
  mockery::expect_args(outcome_1_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(1, 3)))
  mockery::expect_args(outcome_2_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(4)))
})

test_that("hazard resolves three competing outcomes", {
  size <- 4
  population <- individual::Bitset$new(size)$not()

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
    size = size,
    outcomes = list(outcome_1, outcome_2, outcome_3),
    rng = mockery::mock(c(.1, .5, .8, .8))
  )
  
  outcome_1$set_rates(population, c(5, 5, 5, 5))
  outcome_2$set_rates(population, c(5, 5, 5, 5))
  outcome_3$set_rates(population, c(5, 5, 5, 5))
  
  hazard$resolve(0)
   
  mockery::expect_args(outcome_1_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(1)))
  mockery::expect_args(outcome_2_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(2)))
  mockery::expect_args(outcome_3_process, 1, 0,
                       individual::Bitset$new(size)$insert(c(3, 4)))
})
