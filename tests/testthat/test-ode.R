test_that('ODE stays at equilibrium with no biting', {
  foim <- 0
  parameters <- get_parameters()
  parameters$variety_proportions <- 1
  model <- parameterise_ode(parameters, foim)[[1]]
  timesteps <- 365 * 10

  counts <- c()
  
  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, mosquito_model_get_states(model)))
    mosquito_model_step(model, foim)
  }

  expected <- c()
  equilibrium <- initial_mosquito_counts(parameters, foim)
  for (t in seq(timesteps)) {
    expected <- rbind(expected, c(t, equilibrium))
  }

  expect_equal(counts, expected, tolerance=1e-4)
})

test_that('ODE stays at equilibrium with foim > 0', {
  foim <- .1
  parameters <- get_parameters()
  parameters$variety_proportions <- 1
  model <- parameterise_ode(parameters, foim)[[1]]
  timesteps <- 365 * 10

  counts <- c()

  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, mosquito_model_get_states(model)))
    mosquito_model_step(model, foim)
  }

  expected <- c()
  equilibrium <- initial_mosquito_counts(get_parameters(), foim)
  for (t in seq(timesteps)) {
    expected <- rbind(expected, c(t, equilibrium))
  }

  expect_equal(counts, expected, tolerance=1e-4)
})

test_that('Changing FOIM stabilises', {
  foim_0 <- .1
  foim_1 <- .5
  parameters <- get_parameters()
  parameters$variety_proportions <- 1
  model <- parameterise_ode(parameters, foim_0)[[1]]
  timesteps <- 200
  change <- 50

  counts <- c()

  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, mosquito_model_get_states(model)))
    if (t < change) {
      foim <- foim_0
    } else {
      foim <- foim_1
    }
    mosquito_model_step(model, foim)
  }

  initial_eq <- initial_mosquito_counts(parameters, foim_0)
  final_eq <- initial_mosquito_counts(parameters, foim_1)

  expect_equal(counts[1,], c(1, initial_eq), tolerance=1e-4)
  expect_equal(counts[timesteps,], c(timesteps, final_eq), tolerance=1e-4)
  expect_false(any(is.na(counts)))
})
