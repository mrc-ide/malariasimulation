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
   
  mockery::expect_args(
   outcome_1_process,
    1,
    0,
    individual::Bitset$new(size)$insert(c(2, 3))
  )
  mockery::expect_args(
    outcome_2_process,
    1,
    0,
    individual::Bitset$new(size)$insert(c(1, 4))
  )
})

test_that("hazard resolves partial outcomes", {
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
   
  mockery::expect_args(
   outcome_1_process,
    1,
    0,
    individual::Bitset$new(size)$insert(c(3))
  )
  mockery::expect_args(
    outcome_2_process,
    1,
    0,
    individual::Bitset$new(size)$insert(c(2, 4))
  )
})
