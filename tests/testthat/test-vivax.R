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
    c("du",
      "init_ib", "init_iva", "init_ivm", "init_id", # initial immunity parameters
      "rb", "rva", "rid", "rvm", # rates of immune loss
      "ub", "uv", "ud", # immunity non-boosting periods
      "pvm", # maternal immunity parameters
      "b0", "b1", "ib0", "kb", # blood immunity parameters
      "fd0", "ad", "gammad", "d1", "id0", "kd", # asymptomatic detection by light microscopy parameters
      "theta0", "theta1", "iv0", "kv", "fv0", "av", "gammav", # severe immunity parameters
      "gamma1" # asymptomatic infected infectivity towards mosquitos parameter
      )
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

## Test anti-parasite immunity subpatent progression functions
test_that('Test anti-parasite immunity function', {
  vivax_parameters <- get_parameters(parasite = "vivax",
                                     overrides = list(s_proportion = 0,
                                                      d_proportion = 0,
                                                      a_proportion = 0,
                                                      u_proportion = 1,
                                                      t_proportion = 0))
  variables <- create_variables(vivax_parameters)
  index <- variables$state$get_index_of("U")
  dpcr_min <- vivax_parameters$dpcr_min
  dpcr_max <- vivax_parameters$dpcr_max
  apcr50 <- vivax_parameters$apcr50
  kpcr <- vivax_parameters$kpcr
  
  expect_equal(object = anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$iaa$get_values(index),
    variables$iam$get_values(index)),
    expected = rep(dpcr_max,100))
  
  ## Change initial values of ID, and IDM, check they are the same
  variables$iaa <- individual::DoubleVariable$new(1:100)
  UAA_durations <- anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$iaa$get_values(index),
    variables$iam$get_values(index))
  
  variables$iaa <- individual::DoubleVariable$new(rep(0,100))
  variables$iam <- individual::DoubleVariable$new(1:100)
  UAM_durations <- anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$iaa$get_values(index),
    variables$iam$get_values(index))
  
  expect_equal(object = UAA_durations, expected = UAM_durations)
  
  ## Check convergence to min_du at high immunity
  variables$iaa <- individual::DoubleVariable$new(rep(1E7,100))
  variables$iam <- individual::DoubleVariable$new(rep(0,100))
  expect_equal(object = anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$iaa$get_values(index),
    variables$iam$get_values(index)),
    expected = rep(dpcr_min,100),
    tolerance = 1E-2)
})

test_that('that vivax patent prevalence rendering works', {
  
  timestep <- 1
  state <- individual::CategoricalVariable$new(
    c('U', 'A', 'D', 'S', 'Tr'),
    c('U', 'A', 'D', 'D', 'D', 'S')
  )
  birth <- individual::IntegerVariable$new(
    -c(3, 4, 5, 1, 11, 6) * 365
  )
  immunity <- individual::DoubleVariable$new(rep(1, 6))
  vivax_parameters <- get_parameters(parasite = "vivax")
  renderer <- mock_render(1)
  
  process <- create_prevalence_renderer(
    state,
    birth,
    immunity,
    vivax_parameters,
    renderer
  )
  
  process(timestep)
  
  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_detect_lm_730_3650',
    2,
    timestep
  )
  
  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_detect_pcr_730_3650',
    3,
    timestep
  )
  
})
