test_that('1 correlation between rounds gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .2
  parameters <- get_parameters(list(
    human_population = pop,
    pev = TRUE
  ))
  correlations <- get_correlation_parameters(parameters)
  correlations$inter_round_rho('pev', 1)
  round_1 <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  round_2 <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  expect_equal(sum(round_1), pop * .2, tolerance=1e2)
  expect_equal(sum(round_2), pop * .2, tolerance=1e2)
  expect_equal(sum(round_1 & round_2), pop * .2, tolerance=1e2)
})

test_that('0 correlation between rounds gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .5
  parameters <- get_parameters(list(
    human_population = pop,
    pev = TRUE
  ))
  correlations <- get_correlation_parameters(parameters)
  correlations$inter_round_rho('pev', 0)
  round_1 <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  round_2 <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  expect_equal(
    length(intersect(which(round_1), which(round_2))),
    pop * .5,
    tolerance=1e2
  )
  expect_equal(sum(round_1), sum(round_2), tolerance=1e2)
  expect_equal(sum(round_1), pop * .5, tolerance=1e2)
})

test_that('1 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .2
  mda_coverage <- .2
  parameters <- get_parameters(list(
    human_population = pop,
    pev = TRUE,
    mda = TRUE
  ))
  correlations <- get_correlation_parameters(parameters)
  correlations$inter_round_rho('pev', 1)
  correlations$inter_round_rho('mda', 1)
  correlations$inter_intervention_rho('pev', 'mda', 1)
  vaccine_sample <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, correlations)

  expect_equal(sum(vaccine_sample), pop * .2, tolerance=1e2)
  expect_equal(sum(mda_sample), pop * .2, tolerance=1e2)
  expect_equal(sum(vaccine_sample & mda_sample), pop * .2, tolerance=1e2)
})

test_that('0 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .2
  mda_coverage <- .2
  parameters <- get_parameters(list(
    human_population = pop,
    pev = TRUE,
    mda = TRUE
  ))
  correlations <- get_correlation_parameters(parameters)
  correlations$inter_round_rho('pev', 1)
  correlations$inter_round_rho('mda', 1)
  correlations$inter_intervention_rho('pev', 'mda', 0)
  vaccine_sample <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, correlations)
  expect_equal(
    length(intersect(which(vaccine_sample), which(mda_sample))),
    pop * .5,
    tolerance=1e2
  )
  expect_equal(sum(vaccine_sample), sum(mda_sample), tolerance=1e2)
  expect_equal(sum(vaccine_sample), pop * .5, tolerance=1e2)
})

test_that('-1 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  vaccine_coverage <- .2
  mda_coverage <- .2
  parameters <- get_parameters(list(
    human_population = pop,
    pev = TRUE,
    mda = TRUE
  ))
  correlations <- get_correlation_parameters(parameters)
  correlations$inter_round_rho('pev', 1)
  correlations$inter_round_rho('mda', 1)
  correlations$inter_intervention_rho('pev', 'mda', -1)
  vaccine_sample <- sample_intervention(target, 'pev', vaccine_coverage, correlations)
  mda_sample <- sample_intervention(target, 'mda', mda_coverage, correlations)
  expect_equal(length(intersect(which(vaccine_sample), which(mda_sample))), 0)
  expect_equal(sum(vaccine_sample), .2 * pop, tolerance=1e2)
  expect_equal(sum(mda_sample), .2 * pop, tolerance=1e2)
})
