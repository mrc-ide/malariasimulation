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

test_that('correlation between rounds is preserved when adding interventions', {
  pop <- 1e6
  target <- seq(pop)

  pev_coverage_1 <- .2
  pev_coverage_2 <- .4

  initial <- CorrelationParameters$new(pop, c('pev'))
  initial$inter_round_rho('pev', 1)

  restored <- CorrelationParameters$new(pop, c('pev', 'mda'))
  restored$inter_round_rho('pev', 1)
  restored$inter_round_rho('mda', 1)
  restored$inter_intervention_rho('pev', 'mda', 1)
  restored$restore_state(0, initial$save_state())

  round_1 <- sample_intervention(target, 'pev', pev_coverage_1, initial)
  round_2 <- sample_intervention(target, 'pev', pev_coverage_2, restored)

  expect_equal(sum(round_1), pop * pev_coverage_1, tolerance=.1)
  expect_equal(sum(round_2), pop * pev_coverage_2, tolerance=.1)

  expect_equal(
    sum(round_1 & round_2),
    pop * min(pev_coverage_1, pev_coverage_2),
    tolerance=.1)

  expect_equal(
    sum(round_1 | round_2),
    pop * max(pev_coverage_1, pev_coverage_2),
    tolerance=.1)
})

test_that('1 correlation between interventions gives sensible samples when restored', {
  pop <- 1e6
  target <- seq(pop)

  pev_coverage <- .2
  mda_coverage <- .2

  initial <- CorrelationParameters$new(pop, c('pev'))
  initial$inter_round_rho('pev', 1)

  restored <- CorrelationParameters$new(pop, c('pev', 'mda'))
  restored$inter_round_rho('pev', 1)
  restored$inter_round_rho('mda', 1)
  restored$inter_intervention_rho('pev', 'mda', 1)
  restored$restore_state(0, initial$save_state())

  pev_sample <- sample_intervention(target, 'pev', pev_coverage, initial)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, restored)

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

test_that('0 correlation between interventions gives sensible samples when restored', {
  pop <- 1e6
  target <- seq(pop)

  pev_coverage <- .2
  mda_coverage <- .2

  initial <- CorrelationParameters$new(pop, c('pev'))
  initial$inter_round_rho('pev', 1)

  restored <- CorrelationParameters$new(pop, c('pev', 'mda'))
  restored$inter_round_rho('pev', 1)
  restored$inter_round_rho('mda', 1)
  restored$inter_intervention_rho('pev', 'mda', 0)
  restored$restore_state(0, initial$save_state())

  pev_sample <- sample_intervention(target, 'pev', pev_coverage, initial)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, restored)

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

test_that('-1 correlation between interventions gives sensible samples when restored', {
  pop <- 1e6
  target <- seq(pop)
  pev_coverage <- .2
  mda_coverage <- .2

  initial <- CorrelationParameters$new(pop, c('pev'))
  initial$inter_round_rho('pev', 1)

  restored <- CorrelationParameters$new(pop, c('pev', 'mda'))
  restored$inter_round_rho('pev', 1)
  restored$inter_round_rho('mda', 1)
  restored$inter_intervention_rho('pev', 'mda', -1)
  restored$restore_state(0, initial$save_state())

  pev_sample <- sample_intervention(target, 'pev', pev_coverage, initial)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, restored)

  expect_equal(sum(pev_sample), pop * pev_coverage, tolerance=.1)
  expect_equal(sum(mda_sample), pop * mda_coverage, tolerance=.1)

  expect_equal(sum(pev_sample & mda_sample), 0, tolerance=.1)
  expect_equal(
    sum(pev_sample | mda_sample),
    pop * (pev_coverage + mda_coverage),
    tolerance=.1)
})

test_that("rcondmvnorm has correct distribution", {
  set.seed(123)

  pop <- 1e6

  # These are completely arbitrary values. The statistics from the simulated
  # sample will get compared back to this.
  means <- c(-5, 7, 0, 0.3)
  variance <- c(0.3, 0.6, 0.9, 0.2)
  correlation <- matrix(0, ncol=4, nrow=4)
  correlation[1,2] <- correlation[2,1] <- -0.6
  correlation[1,3] <- correlation[3,1] <- 0.4
  correlation[1,4] <- correlation[4,1] <- 1
  correlation[2,3] <- correlation[3,2] <- 0
  correlation[2,4] <- correlation[4,2] <- 0.1
  correlation[3,4] <- correlation[4,3] <- 0.5

  covariance <- outer(variance, variance) * correlation
  diag(covariance) <- variance

  # These are indices of variables that get simulated together.
  wy.ind <- c(1, 3)
  xz.ind <- c(2, 4)

  # Simulate from the MVN for a subset of the dimensions. We intentionally pass
  # in a subset of the mean and covariance matrices, as the rest of the
  # parameters are not needed and may not be known. `dependent.ind` is relative
  # to the subsetted mean and covariance, and therefore is set the first two
  # indices as opposed to `wy.ind`.
  wy <- rcondmvnorm(pop, means[wy.ind], covariance[wy.ind, wy.ind], given=NULL,
                    dependent.ind=c(1,2), given.ind=NULL)
  expect_equal(dim(wy), c(pop, 2))
  expect_equal(apply(wy, 2, mean), means[wy.ind], tolerance=0.001)
  expect_equal(cov(wy), covariance[wy.ind, wy.ind], tolerance=0.001)

  # Now simulate some more, but conditional on the existing values.
  # The call to rcondmvnorm needs all the mean and covariance parameters,
  # including the ones that have already been simulated.
  xz <- rcondmvnorm(pop, means, covariance, given=wy,
                    dependent.ind=xz.ind, given.ind=wy.ind)
  expect_equal(dim(xz), c(pop, 2))
  expect_equal(apply(xz, 2, mean), means[xz.ind], tolerance=0.001)
  expect_equal(cov(xz), covariance[xz.ind, xz.ind], tolerance=0.001)

  # Stitch the variables back together and make sure the covariance across the
  # separately simulated values match the expected one.
  values <- cbind(wy[,1], xz[,1], wy[,2], xz[,2])
  expect_equal(apply(values, 2, mean), means, tolerance=0.001)
  expect_equal(cov(values), covariance, tolerance=0.001)
})

# This would usually be considered an uninteresting edge case, since null
# variance just yields a constant vector. The way we use the function in the
# simulation though, a null variance is actually the default and most common
# case, as it is used where the different rounds of an intervention have no
# correlation.
#
# We need to make sure the other variables are still correctly simulated.
test_that("rcondmvnorm allows null variance of given variables", {
  set.seed(123)

  pop <- 1e6

  means <- c(-5, 7, 0)
  variance <- c(0.3, 0, 0.9)
  correlation <- matrix(0, ncol=3, nrow=3)
  correlation[1,3] <- 0.4
  correlation[3,1] <- 0.4

  covariance <- outer(variance, variance) * correlation
  diag(covariance) <- variance

  # These are indices of variables that get simulated together.
  y.ind <- 2
  xz.ind <- c(1,3)

  # Simulate one dimension. Because its variance was null, the result is a
  # constant vector.
  y <- rcondmvnorm(pop, means, covariance, given=NULL,
                   dependent.ind=y.ind, given.ind=NULL)
  expect_equal(as.vector(y), rep(means[y.ind], pop))

  # Simulate the rest. The original dimension has no influence on the result.
  xz <- rcondmvnorm(pop, means, covariance, given=y,
                    dependent.ind=xz.ind, given.ind=y.ind)

  xyz <- cbind(xz[,1], y, xz[,2])
  expect_equal(apply(xyz, 2, mean), means, tolerance=0.001)
  expect_equal(cov(xyz), covariance, tolerance=0.001)
})
