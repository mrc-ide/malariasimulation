test_that('1 correlation between rounds gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)

  coverage_1 <- .2
  coverage_2 <- .4

  correlations <- CorrelationParameters$new(pop, c('pev'))
  correlations$inter_round_rho('pev', 1)

  round_1 <- sample_intervention(target, 'pev', coverage_1, correlations)
  round_2 <- sample_intervention(target, 'pev', coverage_2, correlations)

  expect_equal(sum(round_1), pop * coverage_1, tolerance=.1)
  expect_equal(sum(round_2), pop * coverage_2, tolerance=.1)

  expect_equal(
    sum(round_1 & round_2),
    pop * min(coverage_1, coverage_2),
    tolerance=.1)

  expect_equal(
    sum(round_1 | round_2),
    pop * max(coverage_1, coverage_2),
    tolerance=.1)
})

test_that('0 correlation between rounds gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)

  coverage_1 <- .2
  coverage_2 <- .4

  correlations <- CorrelationParameters$new(pop, c('pev'))
  correlations$inter_round_rho('pev', 0)

  round_1 <- sample_intervention(target, 'pev', coverage_1, correlations)
  round_2 <- sample_intervention(target, 'pev', coverage_2, correlations)

  expect_equal(sum(round_1), pop * coverage_1, tolerance=.1)
  expect_equal(sum(round_2), pop * coverage_2, tolerance=.1)

  expect_equal(
    sum(round_1 & round_2),
    pop * coverage_1 * coverage_2,
    tolerance=.1)

  expect_equal(
    sum(round_1 | round_2),
    pop * (coverage_1 + coverage_2 - (coverage_1 * coverage_2)),
    tolerance=.1)
})

test_that('1 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)

  pev_coverage <- .2
  mda_coverage <- .4

  correlations <- CorrelationParameters$new(pop, c('pev', 'mda'))
  correlations$inter_round_rho('pev', 1)
  correlations$inter_round_rho('mda', 1)
  correlations$inter_intervention_rho('pev', 'mda', 1)

  pev_sample <- sample_intervention(target, 'pev', pev_coverage, correlations)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, correlations)

  expect_equal(sum(pev_sample), pop * pev_coverage, tolerance=.1)
  expect_equal(sum(mda_sample), pop * mda_coverage, tolerance=.1)

  expect_equal(
    sum(pev_sample & mda_sample),
    pop * min(pev_coverage, mda_coverage),
    tolerance=.1)

  expect_equal(
    sum(pev_sample | mda_sample),
    pop * max(pev_coverage, mda_coverage),
    tolerance=.1)
})

test_that('0 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)

  pev_coverage <- .2
  mda_coverage <- .4

  correlations <- CorrelationParameters$new(pop, c('pev', 'mda'))
  correlations$inter_round_rho('pev', 1)
  correlations$inter_round_rho('mda', 1)
  correlations$inter_intervention_rho('pev', 'mda', 0)

  pev_sample <- sample_intervention(target, 'pev', pev_coverage, correlations)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, correlations)

  expect_equal(sum(pev_sample), pop * pev_coverage, tolerance=.1)
  expect_equal(sum(mda_sample), pop * mda_coverage, tolerance=.1)

  expect_equal(
    sum(pev_sample & mda_sample),
    pop * pev_coverage * mda_coverage,
    tolerance=.1)

  expect_equal(
    sum(pev_sample | mda_sample),
    pop * (pev_coverage + mda_coverage - (pev_coverage * mda_coverage)),
    tolerance=.1)
})

test_that('-1 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)

  pev_coverage <- .2
  mda_coverage <- .4

  correlations <- CorrelationParameters$new(pop, c('pev', 'mda'))
  correlations$inter_round_rho('pev', 1)
  correlations$inter_round_rho('mda', 1)
  correlations$inter_intervention_rho('pev', 'mda', -1)

  pev_sample <- sample_intervention(target, 'pev', pev_coverage, correlations)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, correlations)

  expect_equal(sum(pev_sample), pop * pev_coverage, tolerance=.1)
  expect_equal(sum(mda_sample), pop * mda_coverage, tolerance=.1)

  expect_equal(sum(pev_sample & mda_sample), 0, tolerance=.1)
  expect_equal(
    sum(pev_sample | mda_sample),
    pop * (pev_coverage + mda_coverage),
    tolerance=.1)
})
