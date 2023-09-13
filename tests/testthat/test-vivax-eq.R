test_that('Initial states are consistent with equilibrium', {
  skip_on_ci()
  population <- 100000

  EIRs <- c(1, 5, 10, 100)

  eq_params <- get_parameters(parasite = "vivax", overrides = list(human_population = population))

  expected_states <- sapply(
    EIRs,
    function(EIR) {

      eq <- malariaEquilibriumVivax::human_equilibrium_vivax_full_het(
        EIR = EIR,
        ft = sum(get_treatment_coverages(eq_params, 1)),
        p = eq_params,
        age = EQUILIBRIUM_AGES,
        h = malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)
      )

      colSums(eq$ret[[1]][,c("S","D","A","U","T")] +
        eq$ret[[2]][,c("S","D","A","U","T")]+
        eq$ret[[3]][,c("S","D","A","U","T")]+
        eq$ret[[4]][,c("S","D","A","U","T")]+
        eq$ret[[5]][,c("S","D","A","U","T")])

      eq_summary <- malariaEquilibriumVivax::human_equilibrium_vivax_summarise(eq, eq_params)
      return(colSums(eq_summary$states[,c("S","D","A","U","T")]))
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

  expected_states
  actual_states



})

test_that('Initial immunities are consistent with equilibrium', {
  skip_on_ci()
  population <- 100000

  EIRs <- c(1, 5, 10, 100)

  eq_params <- get_parameters(parasite = "vivax", overrides = list(human_population = population))

  expected_averages <- sapply(
    EIRs,
    function(EIR) {

      eq <- malariaEquilibriumVivax::human_equilibrium_vivax_full_het(
        EIR = EIR,
        ft = sum(get_treatment_coverages(eq_params, 1)),
        p = eq_params,
        age = EQUILIBRIUM_AGES,
        h = malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)
      )

      sum(eq[[1]][,c("ICM")]*eq[[1]][,c("prop")]*malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)$weights[[1]]+
                eq[[2]][,c("ICM")]*eq[[1]][,c("prop")]*malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)$weights[[2]]+
                eq[[3]][,c("ICM")]*eq[[1]][,c("prop")]*malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)$weights[[3]]+
                eq[[4]][,c("ICM")]*eq[[1]][,c("prop")]*malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)$weights[[4]]+
                eq[[5]][,c("ICM")]*eq[[1]][,c("prop")]*malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)$weights[[5]])

      return(colSums(eq$ret[[1]][,c("ID","ICA","ICM","ID","IDM")]*eq$ret[[1]][,c("prop")]*eq$w_het[[1]]+
        eq$ret[[2]][,c("HH","ICA","ICM","ID","IDM")]*eq$ret[[1]][,c("prop")]*eq$w_het[[2]]+
        eq$ret[[3]][,c("HH","ICA","ICM","ID","IDM")]*eq$ret[[1]][,c("prop")]*eq$w_het[[3]]+
        eq$ret[[4]][,c("HH","ICA","ICM","ID","IDM")]*eq$ret[[1]][,c("prop")]*eq$w_het[[4]]+
        eq$ret[[5]][,c("HH","ICA","ICM","ID","IDM")]*eq$ret[[1]][,c("prop")]*eq$w_het[[5]]))

      # eq_summary <- malariaEquilibriumVivax::human_equilibrium_vivax_summarise(eq, eq_params)
      # browser()

      # hh_imm <- Reduce("+",
                       # lapply(1:length(eq$w_het),function(x){eq$ret[[x]][,c("HH","ICA","ICM","ID","IDM","phi","prop")] * eq$w_het[x]}))
      ## Am I weighting prop here?


      # return(colSums(t(t(eq_summary$states[,c("HH","ICA","ICM","ID","IDM")])*eq_summary$states[,c("prop")])))
    })


  actual_averages <- sapply(
    EIRs,
    function(EIR) {
      # browser()
      parms <- set_equilibrium(eq_params, init_EIR = EIR)
      vars <- create_variables(parameters = parms)

      return(c(mean(vars$hypnozoites$get_values()),
               mean(vars$ica$get_values()),
               mean(vars$icm$get_values()),
               mean(vars$id$get_values()),
               mean(vars$idm$get_values())))})

  expected_averages
  actual_averages


})

