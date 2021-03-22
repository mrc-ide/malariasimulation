test_that('ODE stays at equilibrium with a constant total_M', {
  parameters <- get_parameters(list(
    species_proportions = 1
  ))
  total_M <- 1000
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  timesteps <- 365 * 10

  counts <- c()
  
  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, solver_get_states(solvers[[1]])))
    mosquito_model_update(models[[1]], total_M)
    solver_step(solvers[[1]])
  }

  expected <- c()
  equilibrium <- initial_mosquito_counts(parameters)[seq(3)]
  for (t in seq(timesteps)) {
    expected <- rbind(expected, c(t, equilibrium))
  }

  expect_equal(counts, expected, tolerance=1e-4)
})

test_that('Adult ODE stays at equilibrium with a constant foim and mu', {
  foim = 0.5
  parameters <- get_parameters(list(
    species_proportions = 1,
    hybrid_mosquitoes = FALSE,
    init_foim = foim
  ))
  total_M <- 1000
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  timesteps <- 365 * 10

  counts <- c()
  
  for (t in seq(timesteps)) {
    states <- solver_get_states(solvers[[1]])
    counts <- rbind(counts, c(t, states))
    adult_mosquito_model_update(
      models[[1]],
      parameters$mum,
      foim,
      states[ADULT_ODE_INDICES['Pm']]
    )
    solver_step(solvers[[1]])
  }

  expected <- c()
  equilibrium <- initial_mosquito_counts(parameters, parameters$init_foim, total_M)

  for (t in seq(timesteps)) {
    expected <- rbind(expected, c(t, equilibrium))
  }

  #png("ode_S.png")
  #plot(seq(timesteps), counts[, ADULT_ODE_INDICES['Sm'] + 1], type='l', col='red')
  #lines(seq(timesteps), expected[, ADULT_ODE_INDICES['Sm'] + 1], type='l', col='blue')
  #dev.off()

  #png("ode_E.png")
  #plot(seq(timesteps), counts[, ADULT_ODE_INDICES['Pm'] + 1], type='l', col='red')
  #lines(seq(timesteps), expected[, ADULT_ODE_INDICES['Pm'] + 1], type='l', col='blue')
  #dev.off()
  
  #png("ode_I.png")
  #plot(seq(timesteps), counts[, ADULT_ODE_INDICES['Im'] + 1], type='l', col='red')
  #lines(seq(timesteps), expected[, ADULT_ODE_INDICES['Im'] + 1], type='l', col='blue')
  #dev.off()
  
  #png("ode_total_M.png")
  #plot(seq(timesteps), rowSums(counts[, ADULT_ODE_INDICES + 1]), type='l', col='red')
  #lines(seq(timesteps), rowSums(expected[, ADULT_ODE_INDICES + 1]), type='l', col='blue')
  #dev.off()

  expect_equal(counts, expected, tolerance=1e-4)
})

test_that('ODE stays at equilibrium with low total_M', {
  total_M <- 10
  parameters <- get_parameters(list(
    total_M = total_M,
    species = 'all',
    species_proportions = 1
  ))
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  timesteps <- 365 * 10

  counts <- c()

  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, solver_get_states(solvers[[1]])))
    mosquito_model_update(models[[1]], total_M)
    solver_step(solvers[[1]])
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
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  timesteps <- 365 * 10
  change <- 50
  burn_in <- 365 * 5

  counts <- c()

  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, solver_get_states(solvers[[1]])))
    if (t < change) {
      total_M <- total_M_0
    } else {
      total_M <- total_M_1
    }
    mosquito_model_update(models[[1]], total_M)
    solver_step(solvers[[1]])
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
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  timesteps <- 365 * 10

  counts <- c()
  
  for (t in seq(timesteps)) {
    counts <- rbind(counts, c(t, solver_get_states(solvers[[2]])))
    mosquito_model_update(models[[2]], total_M * parameters$species_proportions[[2]])
    solver_step(solvers[[2]])
  }

  expected <- cbind(
    seq(timesteps),
    matrix(0, nrow=timesteps, ncol=3)
  )

  expect_equal(counts, expected)
})
