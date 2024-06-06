#' Test simulation saving and restoring for a given parameter set.
#'
#' This function runs the simulation twice. A first, continuous and uninterrupted
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

  # This function is only used with null correlations. However a null
  # correlation involves sampling random numbers during initialization, which
  # disrupts the global RNG and affects the reproducibility if the size of the
  # matrix is not always the same.
  #
  # We use a single correlation object, that we initialize eagerly, such that
  # the simulation can run undisturbed.
  correlations <- get_correlation_parameters(parameters)
  correlations$mvnorm()

  set.seed(123)
  uninterrupted_run <- run_simulation(
    timesteps,
    parameters = parameters,
    correlations = correlations)

  set.seed(123)
  first_phase <- run_resumable_simulation(
    warmup_timesteps,
    warmup_parameters,
    correlations = correlations)
  second_phase <- run_resumable_simulation(
    timesteps,
    parameters,
    initial_state = first_phase$state,
    restore_random_state = TRUE)

  expect_equal(nrow(first_phase$data), warmup_timesteps)
  expect_equal(nrow(second_phase$data), timesteps - warmup_timesteps)

  # The order of columns isn't always identical, hence why mapequal needs to be
  # used.
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
  # doses do show up, even though they we scheduled during warmup.
  expect_equal(data[data$n_pev_mass_dose_1 > 0, "timestep"], c(100))
  expect_equal(data[data$n_pev_mass_dose_2 > 0, "timestep"], c(70, 145))
  expect_equal(data[data$n_pev_mass_dose_3 > 0, "timestep"], c(115, 190))
})

test_that("TBV intervention can be added when resuming", {
  set_default_tbv <- function(parameters, timesteps) {
    set_tbv(
      parameters,
      timesteps=timesteps,
      coverage=rep(1, length(timesteps)),
      ages=5:60)
  }

  base <- get_parameters()

  data <- test_resume(
    warmup_parameters = base,
    parameters = base %>% set_default_tbv(100))
  expect_equal(data[!is.na(data$n_vaccinated_tbv), "timestep"], 100)

  data <- test_resume(
    warmup_parameters = base %>% set_default_tbv(25),
    parameters = base %>% set_default_tbv(c(25, 100)))
  expect_equal(data[!is.na(data$n_vaccinated_tbv), "timestep"], 100)
})

test_that("MDA intervention can be added when resuming", {
  set_default_mda <- function(parameters, timesteps) {
    parameters %>% set_drugs(list(SP_AQ_params)) %>% set_mda(
      drug = 1,
      timesteps = timesteps,
      coverages = rep(0.8, length(timesteps)),
      min_ages = rep(0, length(timesteps)),
      max_ages = rep(60*365, length(timesteps)))
  }

  base <- get_parameters()

  data <- test_resume(
    warmup_parameters = base,
    parameters = base %>% set_default_mda(100))
  expect_equal(data[data$n_mda_treated > 0, "timestep"], 100)

  data <- test_resume(
    warmup_parameters = base %>% set_default_mda(25),
    parameters = base %>% set_default_mda(c(25, 100)))
  expect_equal(data[data$n_mda_treated > 0, "timestep"], 100)
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

  data <- test_resume(
    warmup_parameters = base,
    parameters = base %>% set_default_bednets(100))
  expect_equal(data[diff(data$n_use_net) > 0, "timestep"], 100)

  data <- test_resume(
    warmup_parameters = base %>% set_default_bednets(25),
    parameters = base %>% set_default_bednets(c(25, 100)))
  expect_equal(data[diff(data$n_use_net) > 0, "timestep"], 100)
})

test_that("Correlations can be set when resuming with new interventions", {
  set.seed(123)

  # When adding a new intervention with a non-zero correlation, we cannot
  # ensure that an uninterrupted run matches the stopped-and-resumed simulation
  # exactly, as the correlation matrix ends up being randomly sampled in a
  # different order. This stops us from using the `test_resume` used throughout
  # the rest of this file. Instead we'll only do stopped-and-resumed simulations
  #' and check its behaviour.
  #
  # We first do a warmup phase with only TBV enabled. We then resume that
  # simulation three times, each time with MDA enabled. Each time we resume the
  # simulation, we set a different correlation parameter between the TBV and
  # MDA interventions, with values -1, 0 and 1.
  #
  # We look at the output data and confirm that the correlation worked as
  # expected. For this we need not only how many people got each intervention,
  # but also how many received both and how many received at least one. This is
  # not normally exposed, so we add an extra process to render these values.
  #
  # For simplicity, for each intervention, we remove any selection process other
  # than overall coverage, such as age range (set to 0-200years) and drug
  # efficacy (set to 100%).
  #
  # We need a large population to make the statistical assertions succeed. We'll
  # only simulate 3 timesteps to keep execution time down: one timestep for
  # warmup during which TBV takes place, one in which MDA takes place and one
  # final timestep to collect the updated variables.
  population <- 10000
  tbv_coverage <- 0.2
  mda_coverage <- 0.4

  warmup_parameters <- get_parameters(overrides=list(human_population=population)) %>%
    set_tbv(
      timesteps=1,
      coverage=tbv_coverage,
      ages=0:200)

  drug <- SP_AQ_params
  drug[1] <- 1. # Override the drug efficacy to 100%
  parameters <- warmup_parameters %>%
    set_drugs(list(drug)) %>%
    set_mda(
      drug = 1,
      timesteps = 2,
      coverages = mda_coverage,
      min_ages = 0,
      max_ages = 200*365)

  create_processes_stub <- function(renderer, variables, events, parameters, ...) {
    p <- function(t) {
      pop <- parameters$human_population
      tbv <- variables$tbv_vaccinated$get_index_of(a=-1, b=0)$not()
      mda <- variables$drug_time$get_index_of(-1)$not()

      renderer$render("total_tbv", tbv$size(), t)
      renderer$render("total_mda", mda$size(), t)
      renderer$render("total_tbv_and_mda", tbv$copy()$and(mda)$size(), t)
      renderer$render("total_tbv_or_mda", tbv$copy()$or(mda)$size(), t)
    }
    c(create_processes(renderer, variables, events, parameters, ...), p)
  }

  mockery::stub(run_resumable_simulation, 'create_processes', create_processes_stub)

  warmup_correlations <- get_correlation_parameters(warmup_parameters)
  warmup_correlations$inter_round_rho('tbv', 1)

  warmup <- run_resumable_simulation(1,
    parameters=warmup_parameters,
    correlations=warmup_correlations)

  zero_correlation <- get_correlation_parameters(parameters)
  zero_correlation$inter_round_rho('tbv', 1)
  zero_correlation$inter_round_rho('mda', 1)

  positive_correlation <- get_correlation_parameters(parameters)
  positive_correlation$inter_round_rho('tbv', 1)
  positive_correlation$inter_round_rho('mda', 1)
  positive_correlation$inter_intervention_rho('tbv', 'mda', 1)

  negative_correlation <- get_correlation_parameters(parameters)
  negative_correlation$inter_round_rho('tbv', 1)
  negative_correlation$inter_round_rho('mda', 1)
  negative_correlation$inter_intervention_rho('tbv', 'mda', -1)

  data <- run_resumable_simulation(
    3,
    initial_state=warmup$state,
    parameters=parameters,
    correlations=zero_correlation)$data %>% tail(1)
  expect_equal(data$total_tbv, population * tbv_coverage, tolerance = 0.1)
  expect_equal(data$total_mda, population * mda_coverage, tolerance = 0.1)
  expect_equal(
    data$total_tbv_and_mda,
    population * (tbv_coverage * mda_coverage),
    tolerance = 0.1)
  expect_equal(
    data$total_tbv_or_mda,
    population * (tbv_coverage + mda_coverage - tbv_coverage * mda_coverage),
    tolerance = 0.1)

  data <- run_resumable_simulation(
    3,
    initial_state=warmup$state,
    parameters=parameters,
    correlations=positive_correlation)$data %>% tail(1)
  expect_equal(data$total_tbv, population * tbv_coverage, tolerance = 0.1)
  expect_equal(data$total_mda, population * mda_coverage, tolerance = 0.1)
  expect_equal(data$total_tbv_and_mda, min(data$total_tbv, data$total_mda))
  expect_equal(data$total_tbv_or_mda, max(data$total_tbv, data$total_mda))

  data <- run_resumable_simulation(
    3,
    initial_state=warmup$state,
    parameters=parameters,
    correlations=negative_correlation)$data %>% tail(1)
  expect_equal(data$total_tbv, population * tbv_coverage, tolerance = 0.1)
  expect_equal(data$total_mda, population * mda_coverage, tolerance = 0.1)
  expect_equal(data$total_tbv_and_mda, 0)
  expect_equal(data$total_tbv_or_mda, data$total_tbv + data$total_mda)
})
