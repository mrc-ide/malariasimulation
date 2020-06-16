test_that('ODE stays at equilibrium with no biting', {
  foim <- 0
  parameters <- get_parameters()
  model <- parameterise_ode(parameters, foim)
  timesteps <- 100

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

  expect_equal(counts, expected)
})

test_that('ODE stays at equilibrium with foim > 0', {
  foim <- .1
  parameters <- get_parameters()
  model <- parameterise_ode(parameters, foim)
  timesteps <- 100

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

  expect_equal(counts, expected)
})
