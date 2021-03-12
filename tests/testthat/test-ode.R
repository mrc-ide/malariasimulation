test_that('ODE stays at equilibrium with a constant total_M', {
  parameters <- get_parameters(list(
    species_proportions = 1
  ))
  total_M <- 1000
  model <- parameterise_ode(parameters)[[1]]
  timesteps <- 365 * 10

  counts <- c()
  
  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, mosquito_model_get_states(model)))
    mosquito_model_step(model, total_M)
  }

  expected <- c()
  equilibrium <- initial_mosquito_counts(parameters)[seq(3)]
  for (t in seq(timesteps)) {
    expected <- rbind(expected, c(t, equilibrium))
  }

  expect_equal(counts, expected, tolerance=1e-4)
})

test_that('ODE stays at equilibrium with low total_M', {
  total_M <- 10
  parameters <- get_parameters(list(
    total_M = total_M,
    species = 'all',
    species_proportions = 1
  ))
  model <- parameterise_ode(parameters)[[1]]
  timesteps <- 365 * 10

  counts <- c()

  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, mosquito_model_get_states(model)))
    mosquito_model_step(model, total_M)
  }

  expected <- c()
  equilibrium <- initial_mosquito_counts(parameters)[seq(3)]
  for (t in seq(timesteps)) {
    expected <- rbind(expected, c(t, equilibrium))
  }

  expect_equal(counts, expected, tolerance=1e-4)
})


test_that('Changing total_M stabilises', {
  total_M_0 <- 500
  total_M_1 <- 400
  parameters <- get_parameters(list(
    total_M = total_M_0,
    species = 'all',
    species_proportions = 1
  ))
  model <- parameterise_ode(parameters)[[1]]
  timesteps <- 365 * 10
  change <- 50
  burn_in <- 365 * 5

  counts <- c()

  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, mosquito_model_get_states(model)))
    if (t < change) {
      total_M <- total_M_0
    } else {
      total_M <- total_M_1
    }
    mosquito_model_step(model, total_M)
  }

  initial_eq <- initial_mosquito_counts(parameters)[seq(3)]
  final_eq <- counts[burn_in, seq(3) + 1]

  expect_equal(counts[1,], c(1, initial_eq), tolerance=1e-4)
  expect_equal(counts[timesteps,], c(timesteps, final_eq), tolerance=1e-4)
  expect_false(isTRUE(all.equal(initial_eq, final_eq)))
  expect_false(any(is.na(counts)))
})

test_that('You can parameterise Total_M = 0 🤯', {
  total_M <- 1000
  parameters <- get_parameters()
  parameters <- set_species(
    parameters,
    species = list(arab_params, fun_params),
    proportions = c(1, 0)
  )
  parameters <- parameterise_total_M(
    parameters,
    total_M
  )
  model <- parameterise_ode(parameters)[[2]]
  timesteps <- 365 * 10

  counts <- c()
  
  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, mosquito_model_get_states(model)))
    mosquito_model_step(model, 0)
  }

  expected <- cbind(
    seq(timesteps),
    matrix(0, nrow=timesteps, ncol=3)
  )

  expect_equal(counts, expected)
})
