test_that('Test age parameter function works', {
  parameters_f <- get_parameters()
  parameters_v <- get_parameters(parasite = "vivax")
  age_limits <- c(0,1,2,3)*365
  parameters_f <- set_epi_outputs(parameters_f,
                                  age_group = age_limits,
                                  incidence = age_limits+1,
                                  clinical_incidence = age_limits+2,
                                  severe_incidence = age_limits+3,
                                  prevalence = age_limits+4
  )
  parameters_v <- set_epi_outputs(parameters_v,
                  age_group = age_limits,
                  incidence = age_limits+1,
                  patent_incidence = age_limits+2,
                  clinical_incidence = age_limits+3,
                  prevalence = age_limits+4,
                  hypnozoite_prevalence = age_limits+5
                  )

  sim_f <- run_simulation(timesteps = 1, parameters_f)
  sim_v <- run_simulation(timesteps = 1, parameters_v)

  expect_true(
    all(
      paste0(rep(c("n_age",
                   "n", "n_inc", "p_inc",
                   "n","n_inc_clinical","p_inc_clinical",
                   "n","n_inc_severe","p_inc_severe",
                   "n","n_detect_pcr","n_detect_lm","p_detect_lm"), each = 3),"_",
             age_limits[-4]+rep(c(0,rep(c(1,2,3), each = 3), rep(4,4)), each = 3),"_",age_limits[-1]-1+rep(c(0,rep(c(1,2,3), each = 3), rep(4,4)), each = 3)
  ) %in%
    names(sim_f)))

  expect_true(all(paste0(rep(c("n_age",
               "n", "n_inc", "p_inc",
               "n", "n_inc_patent", "p_inc_patent",
               "n","n_inc_clinical","p_inc_clinical",
               "n","n_detect_pcr","n_detect_lm",
               "n","n_hypnozoites"), each = 3),"_",
         age_limits[-4]+rep(c(0,rep(c(1,2,3,4), each = 3), rep(5,2)), each = 3),"_",age_limits[-1]-1+rep(c(0,rep(c(1,2,3,4), each = 3), rep(5,2)), each = 3)) %in%
    names(sim_v)))

})
