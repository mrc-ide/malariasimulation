test_that('Initial states are consistent with equilibrium', {
  skip_on_ci()
  population <- 10000
  
  EIRs <- c(1, 10, 100)
  
  eq_params <- get_parameters(parasite = "vivax", overrides = list(human_population = population))
  
  expected_states <- sapply(
    EIRs,
    function(EIR) {
      
      eq <- malariaEquilibriumVivax::vivax_equilibrium(
        age = EQUILIBRIUM_AGES,
        ft = 0,
        EIR = EIR,
        p = translate_vivax_parameters(eq_params))
      
      eq_params <- c(
        list(
          init_foim = eq$FOIM,
          init_EIR = EIR
        ),
        eq_params
      )
      
      eq_params <- parameterise_mosquito_equilibrium(eq_params, EIR)
      
      vars <- create_variables(parameters = eq_params)
      
      eq_expected <- c(sapply(c("S","D","A","U","T"), function(state){sum(eq$states[[state]])}),
                       sapply(c("ICA","ICM","IAA","IAM"), function(x){sum(eq$states$HH * eq$states[[x]])}),
                       sum(colSums(eq$states$HH, dims = 2) * c(1:dim(eq$states$HH)[3]-1)))
      
      observed <- c(c(vars$state$get_size_of("S"),
                      vars$state$get_size_of("D"),
                      vars$state$get_size_of("A"),
                      vars$state$get_size_of("U"),
                      vars$state$get_size_of("Tr"))/population,
                    c(mean(vars$ica$get_values()),
                      mean(vars$icm$get_values()),
                      mean(vars$iaa$get_values()),
                      mean(vars$iam$get_values()),
                      mean(vars$hypnozoites$get_values())))
      
      expect_equal(object = observed, expected =  as.numeric(eq_expected), tolerance = 1E-1)
      
    })
})
