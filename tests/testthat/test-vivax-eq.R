test_that('Initial states are consistent with equilibrium', {
  skip_on_ci()
  population <- 100000

  EIRs <- c(1, 5, 10, 100, 0.1*365)

  eq_params <- get_parameters(parasite = "vivax", overrides = list(human_population = population,
                                                                   sigma_squared = 1.292285))

  expected_states <- sapply(
    EIRs,
    function(EIR) {
      eq <- malariaEquilibriumVivax::vivax_equilibrium(
        age = EQUILIBRIUM_AGES,
        ft = 0,
        EIR = EIR,
        p = translate_vivax_parameters(eq_params))$states

      return(colSums(do.call(rbind, lapply(eq, function(x){colSums(x[,c("S","D","A","U","T")])}))))
    })

  actual_states <- sapply(
    EIRs,
    function(EIR) {
      parms <- set_equilibrium(eq_params, init_EIR = EIR)
      vars <- create_variables(parameters = parms)

      return(c(vars$state$get_size_of("S"),
               vars$state$get_size_of("D"),
               vars$state$get_size_of("A"),
               vars$state$get_size_of("U"),
               vars$state$get_size_of("Tr")))
    })/population

  expect_equal(object = c(expected_states), expected =  c(actual_states), tolerance = 1E-2)

})

test_that('Initial immunities are consistent with equilibrium', {
  skip_on_ci()
  population <- 100000

  EIRs <- c(1, 5, 10, 100, 0.1*365)

  eq_params <- get_parameters(parasite = "vivax", overrides = list(human_population = population,
                                                                   sigma_squared = 1.292285))
  eq_params <- set_species(parameters = eq_params, species = list(kol_params), proportions = 1)

  expected_averages <- sapply(
    EIRs,
    function(EIR) {
      het <- malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)
      eq <- malariaEquilibriumVivax::vivax_equilibrium(
        age = EQUILIBRIUM_AGES,
        ft = 0,
        EIR = EIR,
        p = translate_vivax_parameters(eq_params))$states
      return(colSums(do.call(rbind, lapply(1:length(het$nodes), function(x){colSums(eq[[x]][,c("ICA","ICM","ID","IDM","HH")]*eq[[x]][,"prop"]*het$weights[x])}))))
    })


  actual_averages <- sapply(
    EIRs,
    function(EIR) {

      parms <- set_equilibrium(eq_params, init_EIR = EIR)
      vars <- create_variables(parameters = parms)

      return(c(mean(vars$ica$get_values()),
               mean(vars$icm$get_values()),
               mean(vars$id$get_values()),
               mean(vars$idm$get_values()),
               mean(vars$hypnozoites$get_values())))})

  expect_equal(object = c(expected_states), expected =  c(actual_states), tolerance = 1E-2)

})

test_that('vivax equilibrium works with multiple species', {
  skip_on_ci()

  simparams_pv <- get_parameters(parasite="vivax")
  params_species_pv <- set_species(parameters = simparams_pv,
                                   species = list(arab_params, fun_params, gamb_params),
                                   proportions = c(0.1,0.3,0.6))
  expect_no_error(set_equilibrium(params_species_pv, 6))
})
