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
                    c("du_max","du_min","ku","au50","b","phi0lm","phi1lm","ic0lm","kclm","d0","ca","init_idm","f","gammal")))
})

## Test subpatent progression functions
test_that('Test subpatent duration function works ', {
  vivax_parameters <- get_parameters(parasite = "vivax",
                                     overrides = list(s_proportion = 0,
                                                      d_proportion = 0,
                                                      a_proportion = 0,
                                                      u_proportion = 1,
                                                      t_proportion = 0))
  variables <- create_variables(vivax_parameters)
  index <- variables$state$get_index_of("U")
  du_min <- vivax_parameters$du_min
  du_max <- vivax_parameters$du_max
  au50 <- vivax_parameters$au50
  ku <- vivax_parameters$ku

  expect_equal(object = anti_parasite_immunity(
    du_min, du_max, au50, ku,
    variables$id$get_values(index),
    variables$idm$get_values(index)),
               expected = rep(du_max,100))

  ## Change initial values of ID, and IDM, check they are the same
  variables$id <- individual::DoubleVariable$new(1:100)
  ID_durations <- anti_parasite_immunity(
    du_min, du_max, au50, ku,
    variables$id$get_values(index),
    variables$idm$get_values(index))

  variables$id <- individual::DoubleVariable$new(rep(0,100))
  variables$idm <- individual::DoubleVariable$new(1:100)
  IDM_durations <- anti_parasite_immunity(
    du_min, du_max, au50, ku,
    variables$id$get_values(index),
    variables$idm$get_values(index))

  expect_equal(object = ID_durations, expected = IDM_durations)

  ## Check convergence to min_du at high immunity
  variables$id <- individual::DoubleVariable$new(rep(1E7,100))
  variables$idm <- individual::DoubleVariable$new(rep(0,100))
  expect_equal(object = anti_parasite_immunity(
    du_min, du_max, au50, ku,
    variables$id$get_values(index),
    variables$idm$get_values(index)),
               expected = rep(du_min,100),
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

  # Set all individuals to asymptomatic
  # And ID immunity to not 0 (ID impacts age-specific asymptomatic infectivity)
  parameters <- get_parameters(overrides = list(s_proportion = 0,
                                                d_proportion = 0,
                                                a_proportion = 1,
                                                u_proportion = 0,
                                                t_proportion = 0,
                                                init_id  = 0.5))

  vivax_parameters <- get_parameters(parasite = "vivax",
                                     overrides = list(s_proportion = 0,
                                                      d_proportion = 0,
                                                      a_proportion = 1,
                                                      u_proportion = 0,
                                                      t_proportion = 0,
                                                      init_id  = 0.5))

  # Generate different age structure
  year <- 365
  ages <- round(c(0.083333, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45,
                  50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 200) * year)

  deathrates <- rep(.1, length(ages)) / 365

  dem_parameters <- set_demography(
    parameters,
    agegroups = ages,
    timesteps = 0,
    deathrates = matrix(deathrates, nrow = 1)
  )

  vivax_dem_parameters <- set_demography(
    vivax_parameters,
    agegroups = ages,
    timesteps = 0,
    deathrates = matrix(deathrates, nrow = 1)
  )

  # vivax asymptomatic infectivity should equal ca
  expect_true(all(
    c(create_variables(vivax_parameters)$infectivity$get_values(),
      create_variables(vivax_dem_parameters)$infectivity$get_values())==vivax_parameters$ca))

  # falciparum asymptomatic infectivity should not equal ca
  expect_false(any(
    c(create_variables(parameters)$infectivity$get_values(),
      create_variables(dem_parameters)$infectivity$get_values())==vivax_parameters$ca))

  # falciparum asymptomatic infectivity should change with age structure
  expect_false(any(
    c(create_variables(parameters)$infectivity$get_values() == create_variables(dem_parameters)$infectivity$get_values())))
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
