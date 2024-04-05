test_that('total_M and EIR functions are consistent with equilibrium EIR', {
  skip_on_ci()
  expected_EIR <- c(1, 5, 10, 100)
  population <- 1000
  parameters <- get_parameters(parasite = "vivax",
                               list(
                                 human_population = population,
                                 enable_heterogeneity = FALSE
                               )
  )
  actual_EIR <- vnapply(expected_EIR, function(EIR) {
    parameters <- set_equilibrium(parameters, EIR)

    #set up arguments for EIR calculation
    variables <- create_variables(parameters)
    vector_models <- parameterise_mosquito_models(parameters, 1)
    solvers <- parameterise_solvers(vector_models, parameters)
    estimated_eir <- calculate_eir(1, solvers, variables, parameters, 0)
    mean(estimated_eir) / population * 365
  })

  expect_equal(
    actual_EIR,
    expected_EIR,
    tolerance = 1e-4
  )
})

test_that('total_M and EIR functions are consistent with equilibrium EIR (with het)', {
  skip_on_ci()
  population <- 1000

  expected_EIR <- c(1, 5, 10, 100)

  parameters <- get_parameters(parasite = "vivax",
                               list(human_population = population))

  actual_EIR <- vnapply(expected_EIR, function(EIR) {
    parameters <- set_equilibrium(parameters, EIR)

    variables <- create_variables(parameters)
    vector_models <- parameterise_mosquito_models(parameters, 1)
    solvers <- parameterise_solvers(vector_models, parameters)
    estimated_eir <- calculate_eir(1, solvers, variables, parameters, 0)
    # age <- get_age(variables$birth$get_values(), 0)
    # psi <- unique_biting_rate(age, parameters)
    # omega <- mean(psi)
    mean(estimated_eir) / population * 365
  })

  expect_equal(
    actual_EIR,
    expected_EIR,
    tolerance = 1e-4
  )
})

test_that('FOIM is consistent with equilibrium', {
  skip_on_ci()
  population <- 10000

  EIRs <- c(1, 5, 10, 100)

  eq_params <- get_parameters(parasite = "vivax")

  ages <- 0:800 / 10
  expected_foim <- vnapply(
    EIRs,
    function(EIR) {
      malariaEquilibriumVivax::vivax_equilibrium(
        age = EQUILIBRIUM_AGES,
        ft = 0,
        EIR = EIR,
        p =  translate_vivax_parameters(eq_params))$FOIM
    }
  )

  actual_foim <- vnapply(
    EIRs,
    function(EIR) {
      parameters <- get_parameters(parasite = "vivax", list(human_population = population))
      parameters <- set_equilibrium(parameters, EIR)
      variables <- create_variables(parameters)
      a <- human_blood_meal_rate(1, variables, parameters, 0)
      age <- get_age(variables$birth$get_values(), 0)
      psi <- unique_biting_rate(age, parameters)
      zeta <- variables$zeta$get_values()
      .pi <- human_pi(psi, zeta)
      calculate_foim(a, sum(.pi * variables$infectivity$get_values()), mixing = 1)
    }
  )
  expect_equal(
    expected_foim,
    actual_foim,
    tolerance = 1e-3
  )
})

test_that('phi is consistent with equilibrium at high EIR (no het)', {
  skip_on_ci()
  population <- 1e5

  EIR <- 100
  ages <- 0:800 / 10

  parameters <- get_parameters(parasite = "vivax", list(
    human_population = population,
    enable_heterogeneity = FALSE
  ))
  parameters <- set_equilibrium(parameters, EIR)

  eq <- malariaEquilibriumVivax::vivax_equilibrium(
    age = EQUILIBRIUM_AGES,
    ft = 0,
    EIR = EIR,
    p = translate_vivax_parameters(parameters))$states

  variables <- create_variables(parameters)
  expect_equal(
    mean(
      clinical_immunity(
        variables$ica$get_values(),
        variables$icm$get_values(),
        parameters
      )
    ),
    sum(c(eq$phi_clin * eq$HH)),
    tolerance=1e-2
  )

  # Check phi_patent
  expect_equal(
    mean(
      anti_parasite_immunity(min = parameters$philm_min, max = parameters$philm_max,
                             a50 = parameters$alm50, k = parameters$klm,
                             id = variables$id$get_values(),
                             idm = variables$idm$get_values()
      )
    ),
    sum(c(eq$phi_patent * eq$HH)),
    tolerance=1e-2
  )

  # Check r_PCR
  expect_equal(
    mean(
      1/anti_parasite_immunity(min = parameters$dpcr_min, max = parameters$dpcr_max,
                               a50 = parameters$apcr50, k = parameters$kpcr,
                               id = variables$id$get_values(),
                               idm = variables$idm$get_values()
      )
    ),
    sum(c(eq$r_PCR * eq$HH)),
    tolerance=1e-2
  )

})

test_that('phi is consistent with equilibrium at high EIR', {
  skip_on_ci()
  population <- 1e5

  EIR <- 100
  ages <- 0:999 / 10

  parameters <- get_parameters(parasite = "vivax",
                               list(human_population = population))
  parameters <- set_equilibrium(parameters, EIR)

  eq <- malariaEquilibriumVivax::vivax_equilibrium(
    age = EQUILIBRIUM_AGES,
    ft = 0,
    EIR = EIR,
    p = translate_vivax_parameters(parameters))$states

  variables <- create_variables(parameters)
  het <- statmod::gauss.quad.prob(parameters$n_heterogeneity_groups, dist = "normal")
  expect_equal(
    mean(
      clinical_immunity(
        variables$ica$get_values(),
        variables$icm$get_values(),
        parameters
      )
    ),
    sum(c(eq$phi_clin * eq$HH)),
    tolerance=1e-2
  )

  # Check phi_patent
  expect_equal(
    mean(
      anti_parasite_immunity(min = parameters$philm_min, max = parameters$philm_max,
                             a50 = parameters$alm50, k = parameters$klm,
                             id = variables$id$get_values(),
                             idm = variables$idm$get_values()
      )
    ),
    sum(c(eq$phi_patent * eq$HH)),
    tolerance=1e-2
  )

  # Check r_PCR
  expect_equal(
    mean(
      1/anti_parasite_immunity(min = parameters$dpcr_min, max = parameters$dpcr_max,
                               a50 = parameters$apcr50, k = parameters$kpcr,
                               id = variables$id$get_values(),
                               idm = variables$idm$get_values()
      )
    ),
    sum(c(eq$r_PCR * eq$HH)),
    tolerance=1e-2
  )

})
