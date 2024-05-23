test_that('Test age parameter function works', {
  parameters <- get_parameters()
  age_limits <- c(0,1,2,3)*365
  parameters <- set_epi_outputs(parameters,
                                age_group = age_limits,
                                incidence = age_limits+1,
                                clinical_incidence = age_limits+2,
                                severe_incidence = age_limits+3,
                                prevalence = age_limits+4,
                                ica = age_limits+5,
                                icm = age_limits+6,
                                id = age_limits+7,
                                ib = age_limits+8,
                                iva = age_limits+9,
                                ivm = age_limits+10
  )
  
  sim <- run_simulation(timesteps = 1, parameters)
  
  expect_true(
    all(
      paste0(rep(c("n_age",
                   "n", "n_inc", "p_inc",
                   "n","n_inc_clinical","p_inc_clinical",
                   "n","n_inc_severe","p_inc_severe",
                   "n","n_detect","p_detect",
                   "ica_mean", "icm_mean","id_mean","ib_mean","iva_mean","ivm_mean"), each = 3),"_",
             age_limits[-4]+rep(c(0,rep(c(1:4), each = 3),5:10), each = 3),"_",age_limits[-1]-1+rep(c(0,rep(c(1:4), each = 3),5:10), each = 3)
      ) %in%
        names(sim)))
  
})
