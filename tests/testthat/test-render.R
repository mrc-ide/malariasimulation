test_that('that default rendering works', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  api <- mock_api(
    list(
      human = list(
        U = c(1),
        A = c(2),
        D = c(3),
        S = c(4),
        birth = -c(2, 5, 10, 11) * 365
      )
    ),
    parameters = parameters,
    timestep = 0
  )
  renderer <- create_prevelance_renderer(
    individuals$human,
    states$D,
    states$A,
    variables$birth,
    variables$is_severe
  )
  renderer(api)

  mockery::expect_args(
    api$render,
    1,
    'pv_730_3650',
    2/3
  )
})

test_that('that severe rendering works', {
  year <- 365
  parameters <- get_parameters(list(
    severe_prevalence_rendering_min_ages = c(0, 2) * year,
    severe_prevalence_rendering_max_ages = c(5, 10) * year,
    prevalence_rendering_min_ages = NULL,
    prevalence_rendering_max_ages = NULL
  ))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  api <- mock_api(
    list(
      human = list(
        U = c(1),
        D = c(2, 3),
        S = c(4),
        A = numeric(0),
        birth = -c(2, 5, 10, 11) * 365,
        is_severe = c(0, 1, 0, 0)
      )
    ),
    parameters = parameters,
    timestep = 0
  )
  renderer <- create_prevelance_renderer(
    individuals$human,
    states$D,
    states$A,
    variables$birth,
    variables$is_severe
  )
  renderer(api)

  mockery::expect_args(
    api$render,
    1,
    'pv_severe_0_1825',
    1/2
  )

  mockery::expect_args(
    api$render,
    2,
    'pv_severe_730_3650',
    1/3
  )
})

test_that('that incidence rendering works', {
  year <- 365
  parameters <- get_parameters(list(
    incidence_rendering_min_ages = c(0, 2) * year,
    incidence_rendering_max_ages = c(5, 10) * year,
    prevalence_rendering_min_ages = NULL,
    prevalence_rendering_max_ages = NULL
  ))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  api <- mock_api(
    list(
      human = list(
        birth = -c(2, 5, 10, 11) * 365,
        is_severe = c(0, 1, 0, 0)
      )
    ),
    parameters = parameters,
    timestep = 0
  )
  renderer <- create_incidence_renderer(
    individuals$human,
    variables$birth,
    variables$is_severe
  )

  renderer(api, c(1, 2, 4))

  mockery::expect_args(
    api$render,
    1,
    'inc_0_1825',
    1
  )

  mockery::expect_args(
    api$render,
    2,
    'inc_730_3650',
    2/3
  )
})

test_that('that severe incidence rendering works', {
  year <- 365
  parameters <- get_parameters(list(
    severe_incidence_rendering_min_ages = c(0, 2) * year,
    severe_incidence_rendering_max_ages = c(5, 10) * year,
    prevalence_rendering_min_ages = NULL,
    prevalence_rendering_max_ages = NULL
  ))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  api <- mock_api(
    list(
      human = list(
        birth = -c(2, 6, 10, 11) * 365,
        is_severe = c(0, 1, 0, 0)
      )
    ),
    parameters = parameters,
    timestep = 0
  )
  renderer <- create_incidence_renderer(
    individuals$human,
    variables$birth,
    variables$is_severe
  )

  renderer(api, c(1, 2, 4))

  mockery::expect_args(
    api$render,
    1,
    'inc_severe_0_1825',
    0
  )

  mockery::expect_args(
    api$render,
    2,
    'inc_severe_730_3650',
    1/3
  )
})
