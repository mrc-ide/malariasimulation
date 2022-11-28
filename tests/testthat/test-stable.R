test_that('simulation can run for the equilibrium', {
  skip_on_ci()
  population <- 1e2
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, 10)
  sim <- run_simulation_until_stable(parameters=parameters, tolerance=10)$pre
  expect_gt(nrow(sim), 365)
  expect_equal(mean(tail(sim$EIR_All, 365)) / population, 10, tolerance=10)
})

test_that('simulation can post stability parameters', {
  skip_on_ci()
  population <- 1e2
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, 10)
  post_parameters <- function(warmup) {
    set_bednets(
      parameters,
      timesteps = warmup + c(5, 50),
      coverages = c(.5, .9),
      retention = 40,
      dn0 = matrix(c(.533, .533), nrow=2, ncol=1),
      rn = matrix(c(.56, .56), nrow=2, ncol=1),
      rnm = matrix(c(.24, .24), nrow=2, ncol=1),
      gamman = c(963.6, 963.6)
    )
  }
  sim <- run_simulation_until_stable(
    parameters=parameters,
    post_parameters=post_parameters,
    tolerance=10,
    post_t=100
  )
  expect_gt(nrow(sim$pre), 365)
  expect_equal(nrow(sim$post), 100)
  expect_true(is.null(sim$pre$n_use_net))
  expect_true(!all(sim$post$n_use_net == 0))
})

test_that('simulation exits early for max_t', {
  skip_on_ci()
  population <- 1e2
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, 100)
  expect_error(run_simulation_until_stable(parameters=parameters, max_t=365))
})

test_that('simulation can run for two equilibrium solutions', {
  skip_on_ci()
  population <- 1e2
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, 10)
  parametersets <- list(parameters, parameters)
  mixing <- diag(nrow = 2)
  p_captured <- 1 - diag(nrow = 2)
  sim <- run_metapop_simulation_until_stable(
    parameters=parametersets,
    mixing_tt=1,
    mixing=list(mixing),
    p_captured_tt=1,
    p_captured=list(p_captured),
    p_success=1,
    tolerance=10
  )$pre
  expect_gt(nrow(sim[[1]]), 365)
  expect_equal(mean(tail(sim[[1]]$EIR_All, 365)) / population, 10, tolerance=10)
  expect_gt(nrow(sim[[2]]), 365)
  expect_equal(mean(tail(sim[[2]]$EIR_All, 365)) / population, 10, tolerance=10)
})

test_that('simulation can post stability metapop parameters', {
  skip_on_ci()
  population <- 1e2
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, 10)
  parametersets <- list(parameters, parameters)
  post_parameters <- function(warmup) {
    with_bednets <- set_bednets(
      parameters,
      timesteps = warmup + c(5, 50),
      coverages = c(.5, .9),
      retention = 40,
      dn0 = matrix(c(.533, .533), nrow=2, ncol=1),
      rn = matrix(c(.56, .56), nrow=2, ncol=1),
      rnm = matrix(c(.24, .24), nrow=2, ncol=1),
      gamman = c(963.6, 963.6)
    )
    list(with_bednets, with_bednets)
  }
  mixing <- diag(nrow = 2)
  p_captured <- 1 - diag(nrow = 2)
  sim <- run_metapop_simulation_until_stable(
    parameters=parametersets,
    mixing_tt=1,
    mixing=list(mixing),
    p_captured_tt=1,
    p_captured=list(p_captured),
    p_success=1,
    tolerance=10,
    post_t=100,
    post_parameters=post_parameters
  )
  expect_gt(nrow(sim$pre[[1]]), 365)
  expect_equal(nrow(sim$post[[1]]), 100)
  expect_true(is.null(sim$pre[[1]]$n_use_net))
  expect_true(!all(sim$post[[1]]$n_use_net == 0))
  expect_gt(nrow(sim$pre[[2]]), 365)
  expect_equal(nrow(sim$post[[2]]), 100)
  expect_true(is.null(sim$pre[[2]]$n_use_net))
  expect_true(!all(sim$post[[2]]$n_use_net == 0))
})

test_that('simulation exits early for max_t', {
  skip_on_ci()
  population <- 1e2
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, 100)
  mixing <- diag(nrow = 2)
  p_captured <- 1 - diag(nrow = 2)
  expect_error(
    run_metapop_simulation_until_stable(
      parameters=list(parameters, parameters),
      mixing_tt=1,
      mixing=list(mixing),
      p_captured_tt=1,
      p_captured=list(p_captured),
      p_success=1,
      max_t=365
    )
  )
})


