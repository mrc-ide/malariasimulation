test_that('1 correlation between rounds gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  ints <- c(
    rtss = 1
  )
  parameters <- get_parameters(list(
    human_population = pop,
    rtss = TRUE,
    rtss_coverage = .2
  ))
  c_param <- get_correlation_parameters(parameters)
  c_param$inter_round_rho('rtss', 1)
  round_1 <- sample_intervention(target, 'rtss', parameters$rtss_coverage, c_param)
  round_2 <- sample_intervention(target, 'rtss', parameters$rtss_coverage, c_param)
  expect_equal(sum(round_1), pop * .2, tolerance=1e2)
  expect_equal(sum(round_2), pop * .2, tolerance=1e2)
  expect_equal(sum(round_1 & round_2), pop * .2, tolerance=1e2)
})

test_that('0 correlation between rounds gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  ints <- c(
    rtss = 1
  )
  parameters <- get_parameters(list(
    human_population = pop,
    rtss = TRUE,
    rtss_coverage = .5
  ))
  c_param <- get_correlation_parameters(parameters)
  c_param$inter_round_rho('rtss', 0)
  round_1 <- sample_intervention(target, 'rtss', parameters$rtss_coverage, c_param)
  round_2 <- sample_intervention(target, 'rtss', parameters$rtss_coverage, c_param)
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
  ints <- c(
    rtss = 1,
    mda = 2
  )
  parameters <- get_parameters(list(
    human_population = pop,
    rtss_coverage = .2,
    mda_coverage = .2,
    rtss = TRUE,
    mda = TRUE
  ))
  c_param <- get_correlation_parameters(parameters)
  c_param$inter_round_rho('rtss', 1)
  c_param$inter_round_rho('mda', 1)
  c_param$inter_intervention_rho('rtss', 'mda', 1)
  rtss_sample <- sample_intervention(target, 'rtss', parameters$rtss_coverage, c_param)
  mda_sample <- sample_intervention(target, 'mda', parameters$mda_coverage, c_param)

  expect_equal(sum(rtss_sample), pop * .2, tolerance=1e2)
  expect_equal(sum(mda_sample), pop * .2, tolerance=1e2)
  expect_equal(sum(rtss_sample & mda_sample), pop * .2, tolerance=1e2)
})

test_that('0 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  ints <- c(
    rtss = 1,
    mda = 2
  )
  parameters <- get_parameters(list(
    human_population = pop,
    rtss_coverage = .2,
    mda_coverage = .2,
    rtss = TRUE,
    mda = TRUE
  ))
  c_param <- get_correlation_parameters(parameters)
  c_param$inter_round_rho('rtss', 1)
  c_param$inter_round_rho('mda', 1)
  c_param$inter_intervention_rho('rtss', 'mda', 0)
  rtss_sample <- sample_intervention(target, 'rtss', parameters$rtss_coverage, c_param)
  mda_sample <- sample_intervention(target, 'mda', parameters$mda_coverage, c_param)
  expect_equal(
    length(intersect(which(rtss_sample), which(mda_sample))),
    pop * .5,
    tolerance=1e2
  )
  expect_equal(sum(rtss_sample), sum(mda_sample), tolerance=1e2)
  expect_equal(sum(rtss_sample), pop * .5, tolerance=1e2)
})

test_that('-1 correlation between interventions gives sensible samples', {
  pop <- 1e6
  target <- seq(pop)
  ints <- c(
    rtss = 1,
    mda = 2
  )
  parameters <- get_parameters(list(
    human_population = pop,
    rtss_coverage = .2,
    mda_coverage = .2,
    rtss = TRUE,
    mda = TRUE
  ))
  c_param <- get_correlation_parameters(parameters)
  c_param$inter_round_rho('rtss', 1)
  c_param$inter_round_rho('mda', 1)
  c_param$inter_intervention_rho('rtss', 'mda', -1)
  rtss_sample <- sample_intervention(target, 'rtss', parameters$rtss_coverage, c_param)
  mda_sample <- sample_intervention(target, 'mda', parameters$mda_coverage, c_param)
  expect_equal(length(intersect(which(rtss_sample), which(mda_sample))), 0)
  expect_equal(sum(rtss_sample), .2 * pop, tolerance=1e2)
  expect_equal(sum(mda_sample), .2 * pop, tolerance=1e2)
})
