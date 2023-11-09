test_that('Initial states are consistent with equilibrium', {
  skip_on_ci()
  population <- 100000

  EIRs <- c(1, 5, 10, 100, 0.1*365)

  eq_params <- get_parameters(parasite = "vivax", overrides = list(human_population = population,
                                                                   sigma_squared = 1.292285))

  expected_states <- sapply(
    EIRs,
    function(EIR) {

      # eq <- malariaEquilibriumVivax::human_equilibrium(
      #   EIR = EIR,
      #   ft = sum(get_treatment_coverages(eq_params, 1)),
      #   p = eq_params,
      #   age = EQUILIBRIUM_AGES,
      #   h = malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)
      # )
      #
      het = malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)
      eq <- vivax_equilibrium_init_create_combined(age = EQUILIBRIUM_AGES, ft = 0,
                                                   EIR = EIR,
                                                   p = eq_params,
                                                   K_max = 10,
                                                   use_mid_ages = T,
                                                   malariasimulationoutput = T)$ret


      # return(colSums(eq[,c("S","D","A","U","T")]))
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

  expected_states
  actual_states

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
      # eq <- malariaEquilibriumVivax::human_equilibrium(
      #   EIR = EIR,
      #   ft = sum(get_treatment_coverages(eq_params, 1)),
      #   p = eq_params,
      #   age = EQUILIBRIUM_AGES,
      #   h = malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)
      # )
      het = malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)
      eq <- vivax_equilibrium_init_create_combined(age = EQUILIBRIUM_AGES, ft = 0,
                                                   EIR = EIR,
                                                   p = eq_params,
                                                   K_max = 10,
                                                   use_mid_ages = T,
                                                   malariasimulationoutput = T)$ret
      # colSums(eq$states[,c("ICA","ICM","ID","IDM","HH")] * eq$states[,"prop"])
      # browser()
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

  expected_averages
  actual_averages

})

### OK. So now we've establsiehd that the initialisation is correct, let's make some plots to show the differences between the equilibrium values and the outputs.

# for each EIR
# Get initial SDAUT prevalence
# Run 100 times for 10*365 timesteps
#
# parms <- get_parameters(parasite = "vivax", overrides = list(human_population = 1000))
# eq_results <- lapply(c(0.1,1,10,100), function(EIR){
#
#   # Get equilibria
#   eq <- malariaEquilibriumVivax::human_equilibrium_vivax_full_het(
#     EIR = EIR,
#     ft = sum(get_treatment_coverages(eq_params, 1)),
#     p = eq_params,
#     age = EQUILIBRIUM_AGES,
#     h = malariaEquilibriumVivax::gq_normal(eq_params$n_heterogeneity_groups)
#   )
#
#   eq_summary <- malariaEquilibriumVivax::human_equilibrium_vivax_summarise(eq, parms)
#
#
#   parms_eq <- set_equilibrium(init_EIR = EIR, parms)
#   SimRes <- run_simulation_with_repetitions(timesteps = 10*365, repetitions = 10, parms_eq)
#
#
#   Timestep1 <- SimRes |> filter(timestep ==1) |>
#     group_by(timestep) |>
#     summarise(across(everything(), list(mean = mean)))
#
#   Timestep_eq <- SimRes_sum <- SimRes |>
#     group_by(timestep) |>
#     summarise(across(everything(), list(mean = mean))) |>
#     tail(50) |>
#     colMeans()
#
#   eq_sol <- c(NA,NA,NA,NA,eq_summary$foim, NA,NA,NA,NA,NA,NA,NA,NA,NA,
#     colSums(eq_summary$states[,c("S","D","A","U","T")])*1000,
#     colSums(eq_summary$states[,c("ICA","ICM","ID","IDM","HH")]*eq_summary$states[,c("prop")]),
#     NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
#
#
#   return(data.frame(init_EIR = EIR, eq_sol = eq_sol, t_1 = unlist(Timestep1), mod_eq = unlist(Timestep_eq)))
#
# })
#
# eq_results_out <- lapply(eq_results, function(x){x |>
#     tibble::rownames_to_column(var = 'output')})
#
# eq_comb <- do.call("rbind", eq_results_out) |>
#   tidyr::pivot_longer(cols = c(eq_sol,t_1,mod_eq), names_to = "time", values_to = "value") |>
#   filter(!output %in% c("repetition_mean","timestep"))
#
# eq_plot <- eq_comb |>
#   tidyr::pivot_wider(names_from = output, values_from = value) |>
#   mutate(prev = (D_count_mean+A_count_mean+U_count_mean)/1000) |>
#   tidyr::pivot_longer(cols = 3:35, names_to = "output", values_to = "value") |>
#   mutate(output = factor(output, levels = c(paste0(c("S","D","A","U","Tr"),"_count_mean"),"prev",
#                                             paste0(c("E","L","P","Sm","Pm","Im"),"_All_count_mean"),
#                                             paste0(c("EIR","FOIM","mu","total_M"),"_All_mean"),
#                                             paste0(c("hypnozoites","ica","icm","id","idm"),"_mean_mean"),
#                                             paste0(c("infectivity","n_730_3650","n_detect_730_3650","p_detect_730_3650",
#                                                      "n_bitten","n_new_bite_infections","n_hypnozoites","n_relapses","n_decayed",
#                                                      "n_infections_not_relapse","n_infections",
#                                                      "natural_deaths"),"_mean")
#                                             )),
#          time = factor(time, levels = c("eq_sol","t_1","mod_eq"))) |>
#   # filter(!is.na(value)) |>
#   ggplot() +
#   geom_point(aes(x = time, y = value, color = as.factor(init_EIR))) +
#   facet_wrap(~output, scales = "free_y")
#   # facet_grid(as.factor(output)~as.factor(init_EIR), scales = "free_y")
#   # theme_minimal()
#
# # ggsave("../Malaria_Vivax_Projects/MalariaVivaxProjects/Analysis/Equilibrium_function_across_EIR.pdf", eq_plot,
# #        height = 8, width = 15, limitsize = F)
# #


test_that('vivax equilibrium works with multiple species', {
  skip_on_ci()

  simparams_pv <- get_parameters(parasite="vivax")
  params_species_pv <- set_species(parameters = simparams_pv,
                                   species = list(arab_params, fun_params, gamb_params),
                                   proportions = c(0.1,0.3,0.6))
  expect_no_error(set_equilibrium(params_species_pv, 6))
})
