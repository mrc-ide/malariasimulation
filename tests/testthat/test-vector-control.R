test_that('prob_bitten defaults to 1 with no protection', {
  parameters <- get_parameters(list(human_population = 4))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        net_time = c(-1, -1, -1, -1),
        spray_time = c(-1, -1, -1, -1)
      )
    ),
    timestep = 100,
    parameters = parameters
  )

  expect_equal(
    prob_bitten(individuals, variables, 1, api, parameters),
    list(
      prob_bitten_survives = rep(1, 4),
      prob_bitten = rep(1, 4),
      prob_repelled = rep(0, 4)
    )
  )
})

test_that('prob_bitten correctly calculates net only probabilities', {
  parameters <- get_parameters(list(
    bednets = TRUE,
    gamman = 25
  ))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        net_time = c(-1, 5, 50, 100),
        spray_time = c(-1, -1, -1, -1)
      )
    ),
    timestep = 100,
    parameters = parameters
  )

  expect_mapequal(
    prob_bitten(individuals, variables, 1, api, parameters),
    list(
      prob_bitten_survives = c(1, 0.78640, 0.78640, 0.02723),
      prob_bitten = c(1, 0.78640, 0.78640, 0.02723),
      prob_repelled = c(0, 0.2136, 0.2136, 0.4984)
    )
  )
})

test_that('prob_bitten correctly calculates spraying only probabilities', {
  parameters <- get_parameters(list(spraying = TRUE, gammas = 25))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        net_time = c(-1, -1, -1, -1),
        spray_time = c(-1, 5, 50, 100)
      )
    ),
    timestep = 100,
    parameters = parameters
  )

  expect_equal(
    prob_bitten(individuals, variables, 1, api, parameters),
    list(
      prob_bitten_survives = c(1, 1, 1, 0.175112),
      prob_bitten = c(1, 1, 1, 0.806),
      prob_repelled = c(0, 0, 0, 0.194)
    )
  )
})

test_that('prob_bitten correctly combines spraying and net probabilities', {
  parameters <- get_parameters(list(bednets = TRUE, spraying = TRUE, gamman = 25))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        net_time = c(100, 50, 5, -1),
        spray_time = c(-1, 5, 50, 100)
      )
    ),
    timestep = 100,
    parameters = parameters
  )

  expect_equal(
    prob_bitten(individuals, variables, 1, api, parameters),
    list(
      prob_bitten_survives = c(0.027230, 0.786400, 0.786400, 0.175112),
      prob_bitten = c(0.02723, 0.78640, 0.78640, 0.80600),
      prob_repelled = c(0.4984, 0.2136, 0.2136, 0.1940)
    ),
    tolerance=1e-4
  )
})
