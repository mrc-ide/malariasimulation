test_that('Seasonality correctly affects P', {
  parameters <- get_parameters(list(
    model_seasonality = TRUE,
    g0 = 0.285505,
    g = c(-0.325352, -0.0109352, 0.0779865),
    h = c(-0.132815, 0.104675, -0.013919),
    species_proportions = 1
  ))
  total_M <- 1000
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  timesteps <- 365 * 40

  counts <- c()
  
  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, solver_get_states(solvers[[1]])))
    aquatic_mosquito_model_update(
      models[[1]],
      total_M,
      parameters$blood_meal_rates,
      parameters$mum
    )
    solver_step(solvers[[1]])
  }

  burn_in <- 20

  #sample over the second year to check that vector numbers are changing
  offsets <- c(0, 30*4, 30*8)

  sample_vector_counts <- function(year) {
    counts[offsets + (year - 1) * 365, seq(3) + 1]
  }

  first_sample <- sample_vector_counts(burn_in)
  expect_false(zero_range(first_sample))

  #check that vector numbers are periodic after the burn in
  years <- seq(10) + burn_in
  for (year in years) {
    expect_equal(first_sample, sample_vector_counts(year), tolerance=1e-3)
  }
})

test_that('Rainfall is always > 0', {
  g0 <- 1.639449
  g <- c(-2.1978530, 0.6981486, 0.1225242)
  h <- c(-1.4126293, 1.5910602, -0.7558194)
  expect_true(all(vnapply(seq(365), function(t) rainfall(t, g0, g, h)) > 0))
})
