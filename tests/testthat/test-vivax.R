test_that('Test falciparum switch produces same', {
  parameters_def <- get_parameters()
  parameters_fal <- get_parameters(parasite = "falciparum")
  expect_identical(parameters_def, parameters_fal)
})

test_that('Test vivax model runs', {
  vivax_parameters <- get_parameters(parasite = "vivax")
  sim_res <- run_simulation(timesteps = 100, parameters = vivax_parameters)
  expect_equal(nrow(sim_res), 100)
})

test_that('Test parasite = vivax sets parameters$parasite to vivax', {
  vivax_parameters <- get_parameters(parasite = "vivax")
  expect_identical(vivax_parameters$parasite, "vivax")
})

test_that('Test difference between falciparum and vivax parameter lists', {
  falciparum_parameters <- get_parameters(parasite = "falciparum")
  vivax_parameters <- get_parameters(parasite = "vivax")
  
  in_falciparum_not_vivax <- setdiff(names(falciparum_parameters), names(vivax_parameters))
  in_vivax_not_falciparum <- setdiff(names(vivax_parameters), names(falciparum_parameters))
  
  expect_identical(
    in_falciparum_not_vivax,
    c("init_ib", "rb", "ub", "b0", "b1", "ib0", "kb", # blood immunity parameters
    "gamma1") # asymptomatic infected infectivity towards mosquitos parameter
  )
  
  expect_identical(
    in_vivax_not_falciparum,
    c("dpcr_max", "dpcr_min", "kpcr", "apcr50", # human sub-patent state delay
      "init_iaa", "init_iam", "ra", "ua", # antiparasite immunity parameters)
      "b", # probability of infection given an infectious bite
      "philm_min", "philm_max", "klm", "alm50", # probability of light-microscopy detectable infection parameters
      "ca", # light-microscopy detectable infection infectivity towards mosquitos
      "f", "gammal", "init_hyp", "kmax" # hypnozoite parameters
    )
  )
})
