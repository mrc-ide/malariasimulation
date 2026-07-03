# Run `code` with a fixed seed, restoring the global RNG state afterwards. This
# keeps these tests deterministic without perturbing the RNG state seen by later
# test files (some of which rely on the unseeded global RNG).
with_preserved_seed <- function(seed, code) {
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      if (exists(".Random.seed", envir = globalenv())) {
        rm(".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
  }
  set.seed(seed)
  force(code)
}

test_that('set_forced_eir stores the flag, timesteps and values', {
  parameters <- get_parameters()
  parameters <- suppressWarnings(
    set_forced_eir(parameters, timesteps = c(1, 100), eir = c(5, 20))
  )

  expect_true(parameters$force_EIR)
  expect_equal(parameters$force_EIR_timesteps, c(1, 100))
  expect_equal(parameters$force_EIR_values, c(5, 20))
})

test_that('set_forced_eir warns that the EIR is applied directly', {
  expect_warning(
    get_parameters() |> set_forced_eir(timesteps = 1, eir = 5),
    'bypasses the simulated mosquito population'
  )
})

test_that('set_forced_eir warns when the first change point is after timestep 1', {
  expect_warning(
    get_parameters() |> set_forced_eir(timesteps = c(100, 200), eir = c(5, 10)),
    'first forced EIR change point is at timestep 100'
  )
})

test_that('set_forced_eir validates its inputs', {
  parameters <- get_parameters()

  # mismatched lengths
  expect_error(set_forced_eir(parameters, timesteps = c(1, 100), eir = 5))
  # negative EIR
  expect_error(set_forced_eir(parameters, timesteps = c(1, 100), eir = c(5, -1)))
  # non-positive timestep
  expect_error(set_forced_eir(parameters, timesteps = c(0, 100), eir = c(5, 20)))
  # unsorted timesteps
  expect_error(set_forced_eir(parameters, timesteps = c(100, 1), eir = c(5, 20)))
})

test_that('get_forced_eir is a piecewise-constant step function', {
  parameters <- suppressWarnings(
    get_parameters() |> set_forced_eir(timesteps = c(10, 100), eir = c(5, 20))
  )

  expect_equal(get_forced_eir(parameters, 1), 0)    # before first change point
  expect_equal(get_forced_eir(parameters, 10), 5)   # at first change point
  expect_equal(get_forced_eir(parameters, 50), 5)   # between change points
  expect_equal(get_forced_eir(parameters, 100), 20) # at second change point
  expect_equal(get_forced_eir(parameters, 500), 20) # after last change point
})

test_that('forced_species_eir converts per-adult EIR to the native driver', {
  population <- 1000
  parameters <- suppressWarnings(
    get_parameters(list(human_population = population)) |>
      set_forced_eir(timesteps = 1, eir = 10)
  )

  with_preserved_seed(123, {
    psi <- runif(population)
    # single species, proportion 1. The per-adult EIR is converted to a
    # population-average rate using the same scalar set_equilibrium uses.
    omega <- calculate_population_to_adult_EIR_scalar(parameters, 10)
    expected <- (10 / omega) / 365 * population / mean(psi)
    expect_equal(forced_species_eir(parameters, 1, psi, timestep = 50), expected)

    # before the first change point the forced EIR (and driver) is zero
    expect_equal(forced_species_eir(parameters, 1, psi, timestep = 0), 0)
  })
})

test_that('forced_species_eir splits the total EIR by species proportions', {
  population <- 1000
  proportions <- c(0.7, 0.3)
  parameters <- suppressWarnings(
    get_parameters(list(human_population = population)) |>
      set_species(list(gamb_params, arab_params), proportions) |>
      set_forced_eir(timesteps = 1, eir = 10)
  )

  with_preserved_seed(123, {
    psi <- runif(population)
    eir1 <- forced_species_eir(parameters, 1, psi, timestep = 50)
    eir2 <- forced_species_eir(parameters, 2, psi, timestep = 50)
    # each species carries its share of the total EIR
    expect_equal(eir1 / eir2, proportions[1] / proportions[2])
  })
})

test_that('simulate_bites forces EIR and skips mosquito side effects', {
  population <- 4
  timestep <- 5
  parameters <- suppressWarnings(
    get_parameters(list(human_population = population)) |>
      set_forced_eir(timesteps = 1, eir = 10)
  )

  with_preserved_seed(123, {
    renderer <- individual::Render$new(timestep)
    events <- create_events(parameters)
    variables <- create_variables(parameters)
    variables$zeta <- individual::DoubleVariable$new(c(.2, .3, .5, .9))
    variables$infectivity <- individual::DoubleVariable$new(c(.6, 0, .2, .3))
    age <- c(20, 24, 5, 39) * 365

    eqs_update <- mockery::mock()
    sample_mock <- mockery::mock(c(2, 3))
    pois_mock <- mockery::mock(2)

    mockery::stub(simulate_bites, 'rpois', pois_mock)
    mockery::stub(simulate_bites, 'fast_weighted_sample', sample_mock)
    mockery::stub(simulate_bites, 'adult_mosquito_model_update', eqs_update)

    models <- parameterise_mosquito_models(parameters, timestep)
    solvers <- parameterise_solvers(models, parameters)
    lagged_infectivity <- LaggedValue$new(12.5, .001)
    lagged_eir <- list(LaggedValue$new(12, 10))

    bitten <- simulate_bites(
      renderer,
      solvers,
      models,
      variables,
      events,
      age,
      parameters,
      timestep,
      lagged_infectivity,
      lagged_eir
    )

    expect_equal(bitten$bitten_humans$to_vector(), c(2, 3))

    # the mosquito ODE model must not be updated when the EIR is forced
    mockery::expect_called(eqs_update, 0)

    # the expected number of bites is driven by the forced EIR, not the ODE:
    # the Poisson mean is the forced driver scaled by mean(psi)
    psi <- unique_biting_rate(age, parameters)
    expected_driver <- forced_species_eir(parameters, 1, psi, timestep)
    mockery::expect_args(pois_mock, 1, 1, expected_driver * mean(psi))

    # EIR is rendered, mosquito-only quantities are not
    rendered <- renderer$to_dataframe()
    expect_true('EIR_gamb' %in% names(rendered))
    expect_false('FOIM_gamb' %in% names(rendered))
    expect_false('mu_gamb' %in% names(rendered))
  })
})

test_that('create_processes omits mosquito processes when the EIR is forced', {
  parameters <- suppressWarnings(
    get_parameters() |> set_forced_eir(timesteps = 1, eir = 5)
  )

  with_preserved_seed(123, {
    events <- create_events(parameters)
    variables <- create_variables(parameters)
    vector_models <- parameterise_mosquito_models(parameters, 1)
    solvers <- parameterise_solvers(vector_models, parameters)
    renderer <- individual::Render$new(1)

    processes <- create_processes(
      renderer,
      variables,
      events,
      parameters,
      vector_models,
      solvers
    )

    expect_true('biting_process' %in% names(processes))
    expect_false('solver_process' %in% names(processes))
    expect_false('mosquito_state_renderer' %in% names(processes))
    expect_false('vector_count_renderer' %in% names(processes))
  })
})

test_that('create_processes omits mosquito emergence for individual mosquitoes when forced', {
  parameters <- suppressWarnings(
    get_parameters(list(individual_mosquitoes = TRUE)) |>
      set_forced_eir(timesteps = 1, eir = 5)
  )

  with_preserved_seed(123, {
    events <- create_events(parameters)
    variables <- create_variables(parameters)
    vector_models <- parameterise_mosquito_models(parameters, 1)
    solvers <- parameterise_solvers(vector_models, parameters)
    renderer <- individual::Render$new(1)

    processes <- create_processes(
      renderer,
      variables,
      events,
      parameters,
      vector_models,
      solvers
    )

    expect_true('biting_process' %in% names(processes))
    expect_false('mosquito_emergence_process' %in% names(processes))
    expect_false('solver_process' %in% names(processes))
  })
})

test_that('create_processes keeps mosquito processes when the EIR is not forced', {
  parameters <- get_parameters()

  with_preserved_seed(123, {
    events <- create_events(parameters)
    variables <- create_variables(parameters)
    vector_models <- parameterise_mosquito_models(parameters, 1)
    solvers <- parameterise_solvers(vector_models, parameters)
    renderer <- individual::Render$new(1)

    processes <- create_processes(
      renderer,
      variables,
      events,
      parameters,
      vector_models,
      solvers
    )

    expect_true('solver_process' %in% names(processes))
    expect_true('mosquito_state_renderer' %in% names(processes))
    expect_true('vector_count_renderer' %in% names(processes))
  })
})

test_that('a forced-EIR simulation omits mosquito outputs but reports EIR', {
  parameters <- suppressWarnings(
    get_parameters(list(human_population = 1000)) |>
      set_equilibrium(init_EIR = 5) |>
      set_forced_eir(timesteps = 1, eir = 5)
  )

  sim <- with_preserved_seed(123, {
    run_simulation(timesteps = 100, parameters = parameters)
  })

  expect_true('EIR_gamb' %in% names(sim))
  expect_true('n_bitten' %in% names(sim))
  expect_false('FOIM_gamb' %in% names(sim))
  expect_false('mu_gamb' %in% names(sim))
  expect_false('total_M_gamb' %in% names(sim))
  expect_false('Sm_gamb_count' %in% names(sim))

  # infectious bites are delivered
  expect_gt(sum(sim$n_bitten, na.rm = TRUE), 0)
})

test_that('forcing at the equilibrium EIR keeps prevalence near equilibrium', {
  init_EIR <- 20
  timesteps <- 365

  base <- get_parameters(list(human_population = 5000)) |>
    set_equilibrium(init_EIR = init_EIR)
  forced_params <- suppressWarnings(
    get_parameters(list(human_population = 5000)) |>
      set_equilibrium(init_EIR = init_EIR) |>
      set_forced_eir(timesteps = 1, eir = init_EIR)
  )

  reference <- with_preserved_seed(123, {
    run_simulation(timesteps = timesteps, parameters = base)
  })
  forced <- with_preserved_seed(123, {
    run_simulation(timesteps = timesteps, parameters = forced_params)
  })

  ref_pfpr <- mean(reference$n_detect_lm_730_3650 / reference$n_age_730_3650)
  forced_pfpr <- mean(forced$n_detect_lm_730_3650 / forced$n_age_730_3650)

  # the forced run should track the equilibrium run closely
  expect_equal(forced_pfpr, ref_pfpr, tolerance = 0.15)
})

test_that('a step increase in forced EIR increases prevalence', {
  parameters <- suppressWarnings(
    get_parameters(list(human_population = 5000)) |>
      set_equilibrium(init_EIR = 1) |>
      set_forced_eir(timesteps = c(1, 365), eir = c(1, 100))
  )

  sim <- with_preserved_seed(123, {
    run_simulation(timesteps = 365 * 2, parameters = parameters)
  })
  sim$pfpr <- sim$n_detect_lm_730_3650 / sim$n_age_730_3650

  before <- mean(sim$pfpr[300:365])
  after <- mean(sim$pfpr[(365 * 2 - 65):(365 * 2)])

  expect_gt(after, before)
})

test_that('a forced EIR may be combined with interventions and is itself unchanged', {
  # Interventions are allowed alongside a forced EIR. A vector-only intervention
  # (carrying capacity) has no effect because the mosquito model is not stepped,
  # so the forced EIR is identical with and without it.
  forced <- suppressWarnings(
    get_parameters(list(human_population = 1000)) |>
      set_equilibrium(init_EIR = 5) |>
      set_forced_eir(timesteps = 1, eir = 5)
  )
  forced_with_cc <- suppressWarnings(
    forced |>
      set_carrying_capacity(
        timesteps = 50,
        carrying_capacity_scalers = matrix(2, nrow = 1, ncol = 1)
      )
  )

  base_sim <- with_preserved_seed(123, {
    run_simulation(timesteps = 100, parameters = forced)
  })
  cc_sim <- with_preserved_seed(123, {
    run_simulation(timesteps = 100, parameters = forced_with_cc)
  })

  expect_equal(nrow(cc_sim), 100)
  expect_true('EIR_gamb' %in% names(cc_sim))
  # the forced EIR is unaffected by the (inert) vector intervention
  expect_equal(cc_sim$EIR_gamb, base_sim$EIR_gamb)
})

test_that('a forced-EIR run reproduces a matched equilibrium run at low EIR', {
  # At low (unsaturated) EIR, the realised EIR of a forced run should match a
  # standard set_equilibrium run at the same per-adult init_EIR. This locks in
  # the per-adult unit convention.
  init_EIR <- 2
  timesteps <- 365 * 2
  keep <- 366:timesteps

  base <- get_parameters(list(human_population = 10000)) |>
    set_equilibrium(init_EIR = init_EIR)
  forced_params <- suppressWarnings(
    get_parameters(list(human_population = 10000)) |>
      set_equilibrium(init_EIR = init_EIR) |>
      set_forced_eir(timesteps = 1, eir = init_EIR)
  )

  reference <- with_preserved_seed(1, {
    run_simulation(timesteps = timesteps, parameters = base)
  })
  forced <- with_preserved_seed(1, {
    run_simulation(timesteps = timesteps, parameters = forced_params)
  })

  ref_eir <- mean(reference$EIR_gamb[keep])
  forced_eir <- mean(forced$EIR_gamb[keep])

  expect_equal(forced_eir, ref_eir, tolerance = 0.1)
})

test_that('a forced EIR works for the vivax parasite', {
  parameters <- suppressWarnings(
    get_parameters(list(human_population = 1000), parasite = 'vivax') |>
      set_equilibrium(init_EIR = 5) |>
      set_forced_eir(timesteps = 1, eir = 5)
  )

  sim <- with_preserved_seed(123, {
    run_simulation(timesteps = 100, parameters = parameters)
  })

  expect_equal(nrow(sim), 100)
  expect_true('EIR_gamb' %in% names(sim))
  expect_false('FOIM_gamb' %in% names(sim))
  expect_gt(sum(sim$n_bitten, na.rm = TRUE), 0)
})
