test_that('total_M and EIR functions are consistent with equilibrium EIR', {
  EIR <- 5
  jamie_parameters <- malariaEquilibrium::load_parameter_set()
  h_eq <- malariaEquilibrium::human_equilibrium_no_het(
    EIR,
    0,
    jamie_parameters,
    0:99
  )
  foim <- sum(h_eq[,'inf']*h_eq[,'psi'])
  population <- 100000
  parameters <- get_parameters(c(
    translate_jamie(remove_unused_jamie(jamie_parameters)),
    list(
      init_foim = foim,
      variety_proportions = 1,
      human_population = population
    )
  ))
  parameters <- parameterise_mosquito_equilibrium(parameters, EIR)
  m_eq <- initial_mosquito_counts(parameters, foim)

  age <- rep(seq(100) - 1, rowSums(h_eq[,c('S', 'U', 'A', 'D')]) * population)
  xi <- rep(1, length(age))
  infectivity <- m_eq[[6]] * parameters$blood_meal_rate
  expect_equal(
    mean(eir(age, xi, infectivity, parameters)) * 365,
    EIR,
    tolerance = 1
  )
})

test_that('total_M and EIR functions are consistent with equilibrium EIR (with het)', {
  population <- 100000

  EIR <- 50
  n_groups <- 10

  jamie_parameters <- malariaEquilibrium::load_parameter_set()
  het_groups <- statmod::gauss.quad.prob(n_groups, dist = "normal")
  h_eqs <- list()
  foim <- 0
  all_age <- c()
  all_xi <- c()

  for (h in seq(n_groups)) {
    h_eq <- malariaEquilibrium::human_equilibrium_no_het(
      EIR,
      0,
      jamie_parameters,
      0:99
    )
    xi <- exp(
      -jamie_parameters$s2 * .5 + sqrt(
        jamie_parameters$s2
      ) * het_groups$nodes[h]
    )
    age <- rep(
      seq(100) - 1,
      rowSums(h_eq[,c('S', 'U', 'A', 'D')]) * het_groups$weights[h] * population
    )
    all_age <- c(all_age, age)
    all_xi <- c(all_xi, rep(xi, length(age)))
    foim <- foim + sum(h_eq[,'inf'] * h_eq[,'psi']) * het_groups$weights[h] * xi
  }

  parameters <- get_parameters(c(
    translate_jamie(remove_unused_jamie(jamie_parameters)),
    list(
      init_foim = foim,
      variety_proportions = 1,
      human_population = population
    )
  ))
  parameters <- parameterise_mosquito_equilibrium(parameters, EIR)
  m_eq <- initial_mosquito_counts(parameters, foim)
  infectivity <- m_eq[[6]] * parameters$blood_meal_rate

  expect_equal(
    mean(eir(all_age, all_xi, infectivity, parameters)) * 365,
    EIR,
    tolerance = 1
  )
})

test_that('mosquito_limit is set to 0 for 0 EIR', {
  parameters <- parameterise_mosquito_equilibrium(get_parameters(), 0)
  expect_equal(parameters$mosquito_limit, 0)
})

test_that('mosquito_limit is set to a sensible level', {
  EIRs <- c(5, 50, 1000)

  seasonalities <- list(
    list(
      g0 = 7.2,
      g = c(-10, 4.06, -1.08),
      h = c(-4.6, 3.59, -1.32)
    ),
    list(
      g0 = 7.2,
      g = c(-10, 4.06, -1.08),
      h = c(-4.6, 3.59, -1.32)
    )
  )

  for (EIR in EIRs) {
    for (seasonality in seasonalities) {
      parameters <- get_parameters(c(
        seasonality,
        model_seasonality = TRUE,
        init_foim = .1
      ))
      parameters <- parameterise_mosquito_equilibrium(parameters, EIR)
      run_simulation(365, parameters)
    }
  }
  expect_true(TRUE)
})
