test_that('1 correlation between rounds gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .2
  correlations <- CorrelationParameters$new(pop, c('pev'))
  correlations$inter_round_rho('pev', 1)
  round_1 <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  round_2 <- sample_intervention(target, 'pev', vaccine_coverage, correlations)

  expect_equal(sum(round_1), pop * .2, tolerance=.1)
  expect_equal(sum(round_2), pop * .2, tolerance=.1)
  expect_equal(sum(round_1 & round_2), pop * .2, tolerance=.1)
  expect_equal(sum(round_1 | round_2), pop * .2, tolerance=.1)
})

test_that('0 correlation between rounds gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .5
  correlations <- CorrelationParameters$new(pop, c('pev'))
  correlations$inter_round_rho('pev', 0)
  round_1 <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  round_2 <- sample_intervention(target, 'pev', vaccine_coverage, correlations)

  expect_equal(sum(round_1), sum(round_2), tolerance=.1)
  expect_equal(sum(round_1), pop * .5, tolerance=.1)
  expect_equal(sum(round_1 & round_2), pop * .25, tolerance=.1)
  expect_equal(sum(round_1 | round_2), pop * .75, tolerance=.1)
})

test_that('1 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .2
  mda_coverage <- .2
  correlations <- CorrelationParameters$new(pop, c('pev', 'mda'))
  correlations$inter_round_rho('pev', 1)
  correlations$inter_round_rho('mda', 1)
  correlations$inter_intervention_rho('pev', 'mda', 1)
  vaccine_sample <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, correlations)

  expect_equal(sum(vaccine_sample), pop * .2, tolerance=.1)
  expect_equal(sum(mda_sample), pop * .2, tolerance=.1)
  expect_equal(sum(vaccine_sample & mda_sample), pop * .2, tolerance=.1)
  expect_equal(sum(vaccine_sample | mda_sample), pop * .2, tolerance=.1)
})

test_that('0 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .2
  mda_coverage <- .2
  correlations <- CorrelationParameters$new(pop, c('pev', 'mda'))
  correlations$inter_round_rho('pev', 1)
  correlations$inter_round_rho('mda', 1)
  correlations$inter_intervention_rho('pev', 'mda', 0)

  vaccine_sample <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, correlations)

  expect_equal(sum(vaccine_sample), sum(mda_sample), tolerance=.1)
  expect_equal(sum(vaccine_sample), pop * .2, tolerance=.1)
  expect_equal(sum(vaccine_sample & mda_sample), pop * .04, tolerance=.1)
  expect_equal(sum(vaccine_sample | mda_sample), pop * .36, tolerance=.1)
})

test_that('-1 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .2
  mda_coverage <- .2
  correlations <- CorrelationParameters$new(pop, c('pev', 'mda'))
  correlations$inter_round_rho('pev', 1)
  correlations$inter_round_rho('mda', 1)
  correlations$inter_intervention_rho('pev', 'mda', -1)

  vaccine_sample <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, correlations)

  expect_equal(sum(vaccine_sample), .2 * pop, tolerance=.1)
  expect_equal(sum(mda_sample), .2 * pop, tolerance=.1)
  expect_equal(sum(vaccine_sample & mda_sample), 0, tolerance=.1)
  expect_equal(sum(vaccine_sample | mda_sample), pop * .4, tolerance=.1)
})
