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

  expect_true(all(names(falciparum_parameters)[!names(falciparum_parameters) %in% names(vivax_parameters)] %in%
                    c("du","rvm","rva","rb","b0","b1","ib0","kb","theta0","theta1","kv","fv0","av","gammav","iv0","fd0","ad","gammad","d1","id0","kd","ub","uv","gamma1","pvm","init_ivm","init_ib","init_iva")))
  expect_true(all(names(vivax_parameters[!names(vivax_parameters) %in% names(falciparum_parameters)]) %in%
                    c("dpcr_max","dpcr_min","kpcr","apcr50","init_idm","b","philm_min","philm_max","klm","alm50","ca","f","gammal","init_hyp")))
})

## Test subpatent progression functions
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
    variables$id$get_values(index),
    variables$idm$get_values(index)),
    expected = rep(dpcr_max,100))

  ## Change initial values of ID, and IDM, check they are the same
  variables$id <- individual::DoubleVariable$new(1:100)
  ID_durations <- anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$id$get_values(index),
    variables$idm$get_values(index))

  variables$id <- individual::DoubleVariable$new(rep(0,100))
  variables$idm <- individual::DoubleVariable$new(1:100)
  IDM_durations <- anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$id$get_values(index),
    variables$idm$get_values(index))

  expect_equal(object = ID_durations, expected = IDM_durations)

  ## Check convergence to min_du at high immunity
  variables$id <- individual::DoubleVariable$new(rep(1E7,100))
  variables$idm <- individual::DoubleVariable$new(rep(0,100))
  expect_equal(object = anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$id$get_values(index),
    variables$idm$get_values(index)),
    expected = rep(dpcr_min,100),
    tolerance = 1E-2)
})

test_that('Test default vivax incidence rendering works', {

  timestep <- 0
  year <- 365
  birth <- individual::IntegerVariable$new(
    -c(2, 5, 10, 11) * year
  )
  vivax_parameters <- get_parameters(
    parasite = "vivax")

  renderer <- mock_render(1)
  incidence_renderer(
    birth,
    renderer,
    individual::Bitset$new(4)$insert(c(1, 2, 4)),
    individual::Bitset$new(4)$insert(seq(4)),
    c(.1, .2, .3, .4),
    'inc_patent_',
    c(0, 2) * year,
    c(5, 10) * year,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_0_1825',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_inc_patent_0_1825',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    3,
    'p_inc_patent_0_1825',
    0.3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    4,
    'n_730_3650',
    3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    5,
    'n_inc_patent_730_3650',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    6,
    'p_inc_patent_730_3650',
    .6,
    timestep
  )
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

  process <- create_prevelance_renderer(
    state,
    birth,
    immunity,
    vivax_parameters,
    renderer
  )

  mockery::stub(process, 'probability_of_detection', mockery::mock(.5))
  mockery::stub(process, 'bernoulli_multi_p', mockery::mock(1))
  process(timestep)

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_730_3650',
    4,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_detect_pcr_730_3650',
    3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    3,
    'n_detect_lm_730_3650',
    2,
    timestep
  )
})


test_that('Test age structure should not change vivax infectivity', {

  falc_parameters <- get_parameters(
    overrides = list(
      human_population = 1,
      init_id  = 0.5))

  vivax_parameters <- get_parameters(
    parasite = "vivax",
    overrides = list(
      human_population = 1,
      init_id  = 0.5))

  state_mock <- mockery::mock('A', cycle = T)
  mockery::stub(create_variables, 'initial_state', state_mock)

  ages_mock <- mockery::mock(365, cycle = T)
  mockery::stub(create_variables, 'calculate_initial_ages', ages_mock)

  falc_variables <- create_variables(falc_parameters)
  vivax_variables <- create_variables(vivax_parameters)

  expect_equal(falc_variables$infectivity$get_values(), 0.06761596)
  expect_equal(vivax_variables$infectivity$get_values(), 0.1)

  ages_mock <- mockery::mock(365*70, cycle = T)
  mockery::stub(create_variables, 'calculate_initial_ages', ages_mock)

  falc_variables <- create_variables(falc_parameters)
  vivax_variables <- create_variables(vivax_parameters)

  expect_equal(falc_variables$infectivity$get_values(), 0.03785879)
  expect_equal(vivax_variables$infectivity$get_values(), 0.1)

})

test_that('vivax schedule_infections correctly schedules new infections', {
  parameters <- get_parameters(list(human_population = 20), parasite = "vivax")
  variables <- create_variables(parameters)

  infections <- individual::Bitset$new(20)$insert(1:20)
  patent_infections <- individual::Bitset$new(20)$insert(1:16)
  clinical_infections <- individual::Bitset$new(20)$insert(5:15)
  treated <- individual::Bitset$new(20)$insert(7:12)

  infection_mock <- mockery::mock()
  mockery::stub(schedule_infections, 'update_infection', infection_mock)

  schedule_infections(
    variables = variables,
    patent_infections = patent_infections,
    clinical_infections = clinical_infections,
    treated = treated,
    infections = infections,
    parameters = parameters,
    timestep = 42
  )

  actual_infected <- mockery::mock_args(infection_mock)[[1]][[5]]$to_vector()
  actual_asymp_infected <- mockery::mock_args(infection_mock)[[2]][[5]]$to_vector()
  actual_subpatent_infected <- mockery::mock_args(infection_mock)[[3]][[5]]$to_vector()

  expect_equal(
    actual_infected,
    c(5, 6, 13, 14, 15)
  )

  expect_equal(
    actual_asymp_infected,
    c(1, 2, 3, 4, 16)
  )

  expect_equal(
    actual_subpatent_infected,
    c(17, 18, 19, 20)
  )

})

test_that('relapses are recognised', {
  timestep <- 50
  parameters <- get_parameters(parasite = "vivax")

  variables <- list(
    state = individual::CategoricalVariable$new(
      c('D', 'S', 'A', 'U', 'Tr'),
      c('D', 'S', 'A', 'U')
    ),
    drug = individual::DoubleVariable$new(c(0, 0, 0, 0)),
    drug_time = individual::DoubleVariable$new(c(-1, -1, -1, -1)),
    pev_timestep = individual::DoubleVariable$new(c(-1, -1, -1, -1)),
    pev_profile = individual::IntegerVariable$new(c(-1, -1, -1, -1)),
    hypnozoites = individual::IntegerVariable$new(c(0, 1, 2, 3))
  )

  bernoulli_mock <- mockery::mock(3, cycle = TRUE)
  mockery::stub(calculate_infections, 'bernoulli_multi_p', bernoulli_mock)
  bitten_humans <- individual::Bitset$new(4)$insert(c(1, 2, 3, 4))
  renderer <- mock_render(1)

  infections <- calculate_infections(
    variables,
    bitten_humans,
    parameters,
    renderer,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_new_bite_infections',
    1,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_relapses',
    1,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    3,
    'n_infections',
    2,
    timestep
  )
})

test_that('relapses are recognised', {
  timestep <- 50
  parameters <- get_parameters(parasite = "vivax")

  variables <- list(
    state = individual::CategoricalVariable$new(
      c('D', 'S', 'A', 'U', 'Tr'),
      c('D', 'S', 'A', 'U')
    ),
    drug = individual::DoubleVariable$new(c(0, 0, 0, 0)),
    drug_time = individual::DoubleVariable$new(c(-1, -1, -1, -1)),
    pev_timestep = individual::DoubleVariable$new(c(-1, -1, -1, -1)),
    pev_profile = individual::IntegerVariable$new(c(-1, -1, -1, -1)),
    hypnozoites = individual::IntegerVariable$new(c(0, 1, 2, 3))
  )

  bernoulli_mock <- mockery::mock(3, cycle = TRUE)
  mockery::stub(calculate_infections, 'bernoulli_multi_p', bernoulli_mock)
  bitten_humans <- individual::Bitset$new(4)$insert(c(1, 2, 3, 4))
  renderer <- mock_render(1)

  infections <- calculate_infections(
    variables,
    bitten_humans,
    parameters,
    renderer,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_new_bite_infections',
    1,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_relapses',
    1,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    3,
    'n_infections',
    2,
    timestep
  )
})
