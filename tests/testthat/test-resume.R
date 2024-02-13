#' Test simulation saving and restoring for a given parameter set.
#'
#' This function runs the simulation twice. A first, continous and uninterrupted
#' run of the simulation is used as a reference. The second run is split into
#' two phases. Between the two phases, the simulation state is saved and
#' restored. Optionally, the initial warmup phase can use a different set of
#' parameters, by specifying a value for warmup_parameters.
test_resume <- function(
  parameters,
  timesteps = 200,
  warmup_parameters = parameters,
  warmup_timesteps = 50
  ) {
  set.seed(123)
  uninterrupted_run <- run_simulation(timesteps, parameters=parameters)

  set.seed(123)
  first_phase <- run_resumable_simulation(warmup_timesteps, warmup_parameters)
  second_phase <- run_resumable_simulation(
    timesteps,
    parameters,
    initial_state=first_phase$state,
    restore_random_state=TRUE)

  expect_equal(nrow(first_phase$data), warmup_timesteps)
  expect_equal(nrow(second_phase$data), timesteps - warmup_timesteps)

  expect_mapequal(
    second_phase$data,
    uninterrupted_run[-(1:warmup_timesteps),])

  invisible(second_phase$data)
}

test_that('Simulation can be resumed', {
  test_resume(get_parameters(overrides=list(
    individual_mosquitoes = FALSE,
    model_seasonality = TRUE
  )))
  test_resume(get_parameters(overrides=list(
    individual_mosquitoes = TRUE,
    model_seasonality = TRUE
  )))
  test_resume(get_parameters(overrides=list(
    individual_mosquitoes = FALSE,
    model_seasonality = TRUE
  )))
  test_resume(get_parameters(overrides=list(
    individual_mosquitoes = TRUE,
    model_seasonality = TRUE
  )))
})

test_that('PEV intervention can be added when resuming', {
  set_default_mass_pev <- function(parameters, timesteps) {
    n <- length(timesteps)
    set_mass_pev(
      parameters,
      profile = rtss_profile,
      timesteps = timesteps,
      coverages = rep(0.5, n),
      min_wait = 5,
      min_ages = 365*10,
      max_ages = 365*60,
      booster_spacing = NULL,
      booster_coverage = NULL,
      booster_profile = NULL)
  }

  base <- get_parameters(overrides=list(pev_doses=c(0, 45, 90)))
  data <- test_resume(
    warmup_parameters = base,
    parameters = base %>% set_default_mass_pev(100))

  expect_equal(data[data$n_pev_mass_dose_1 > 0, "timestep"], 100)
  expect_equal(data[data$n_pev_mass_dose_2 > 0, "timestep"], 145)
  expect_equal(data[data$n_pev_mass_dose_3 > 0, "timestep"], 190)

  # Add a second mass PEV intervention when resuming the simulation.
  data <- test_resume(
    warmup_parameters = base %>% set_default_mass_pev(25),
    parameters = base %>% set_default_mass_pev(c(25, 100)))

  # The first dose, at time step 25, happens during the warmup and is not
  # returned by test_resume, hence why we don't see it in the data. Follow-up
  # doses do show up, event though they we scheduled during warmup.
  expect_equal(data[data$n_pev_mass_dose_1 > 0, "timestep"], c(100))
  expect_equal(data[data$n_pev_mass_dose_2 > 0, "timestep"], c(70, 145))
  expect_equal(data[data$n_pev_mass_dose_3 > 0, "timestep"], c(115, 190))
})

test_that("TBV intervention can be added when resuming", {
  base <- get_parameters()
  data <- test_resume(
    warmup_parameters = base,
    parameters = base %>% set_tbv(timesteps=100, coverage=1, ages=5:60))
  expect_equal(data[!is.na(data$n_vaccinated_tbv), "timestep"], 100)

  test_resume(
    warmup_parameters = base %>% set_tbv(
      timesteps=25,
      coverage=1,
      ages=5:60),
    parameters = base %>% set_tbv(
      timesteps=c(25, 100),
      coverage=rep(1, 2),
      ages=rep(5:60, 2)))
})

test_that("MDA intervention can be added when resuming", {
  base <- get_parameters()
  data <- test_resume(
    warmup_parameters = base,
    parameters = base %>% set_drugs(list(SP_AQ_params)) %>% set_mda(
      drug = 1,
      timesteps = c(100, 150),
      coverages = rep(0.8, 2),
      min_ages = rep(0, 2),
      max_ages = c(60*365, 2)))

  expect_equal(data[data$n_mda_treated > 0, "timestep"], c(100))
})

test_that("Bednets intervention can be added when resuming", {
  set_default_bednets <- function(parameters, timesteps) {
    n <- length(timesteps)
    set_bednets(
      parameters,
      timesteps = timesteps,
      coverages = rep(0.5, n),
      retention = 25,
      dn0 = matrix(rep(0.533, n), ncol=1),
      rn = matrix(rep(0.56, n), ncol=1),
      rnm = matrix(rep(0.24, n), ncol=1),
      gamman = rep(2.64 * 365, n))
  }

  base <- get_parameters()

  test_resume(
    warmup_parameters = base,
    parameters = base %>% set_default_bednets(100))

  test_resume(
    warmup_parameters = base %>% set_default_bednets(25),
    parameters = base %>% set_default_bednets(c(25, 100)))
})
