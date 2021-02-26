test_that('Seasonality correctly affects P', {
  parameters <- get_parameters(list(
    model_seasonality = TRUE,
    g0 = 0.285505,
    g = c(-0.325352, -0.0109352, 0.0779865),
    h = c(-0.132815, 0.104675, -0.013919),
    species_proportions = 1
  ))
  total_M <- 1000
  model <- parameterise_ode(parameters)[[1]]
  timesteps <- 365 * 40

  counts <- c()
  
  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, mosquito_model_get_states(model)))
    mosquito_model_step(model, total_M)
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
