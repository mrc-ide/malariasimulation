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
  
  prefixes <- c("n_age",
                "n_inc", "p_inc",
                "n_inc_clinical","p_inc_clinical",
                "n_inc_severe","p_inc_severe",
                "n_detect_lm","p_detect_lm","n_detect_pcr",
                "ica_mean", "icm_mean","id_mean","ib_mean","iva_mean","ivm_mean")
  offsets <- c(0, rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 3), 5:10)
  expect_equal(length(prefixes), length(offsets))

  expected <- paste0(rep(prefixes, each = 3),
                     "_", age_limits[-4]+rep(offsets, each = 3),
                     "_", age_limits[-1]-1+rep(offsets, each = 3))
  expect_in(expected, names(sim))
  
  
  expect_error(set_epi_outputs(parameters, iaa = age_limits))
  expect_error(get_parameters(parasite = "vivax") |> 
                 set_epi_outputs(iva = age_limits))
  
})
