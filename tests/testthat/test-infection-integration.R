
test_that('human infection_process creates the correct updates', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events)

  infection_process <- create_infection_process(
    individuals,
    states,
    variables,
    events
  )

  parameters$severe_enabled = 1

  api <- mock_api(
    list(
      human = list(
        S = c(2),
        U = c(3),
        A = c(4),
        D = c(1),
        age = c(20, 24, 5, 39),
        IB = c(.2, .3, .5, .9),
        xi = c(.2, .3, .5, .9),
        ICA = c(.2, .3, .5, .9),
        IVA = c(.2, .3, .5, .9),
        ICM = c(.2, .3, .5, .9),
        IVM = c(.2, .3, .5, .9),
        last_infected = c(-1, -1, 1, -1),
        last_bitten = c(-1, 1, 1, -1),
        ID = c(.2, .3, .5, .9)
      ),
      mosquito = list(
        Im = 1:100,
        variety = c(rep(1, 25), rep(2, 25), rep(3, 50))
      )
    ),
    timestep = 5,
    parameters = parameters
  )

  mockery::stub(
    infection_process,
    'bernoulli',
    mockery::mock(
      c(TRUE, FALSE, TRUE),
      c(TRUE, FALSE),
      c(FALSE, TRUE),
      c(FALSE),
    )
  )

  api$get_scheduled = mockery::mock(4, 1)

  infection_process(api)

  updates <- mockery::mock_args(api$queue_variable_update)

  expect_variable_update(updates[[1]], 'IB', c(.3, 1.9), c(2, 4))
  expect_variable_update(updates[[2]], 'last_bitten', 5, c(2, 4))
  expect_variable_update(updates[[3]], 'ICA', 1.9, 4)
  expect_variable_update(updates[[4]], 'IVA', 1.9, 4)
  expect_variable_update(updates[[5]], 'ID', 1.9, 4)
  expect_variable_update(updates[[6]], 'last_infected', 5, 4)
})

test_that('mosquito_infection_process creates the correct updates', {
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  events <- create_events()
  individuals <- create_individuals(states, variables, create_events())
  mosquito_infection_process <- create_mosquito_infection_process(
    individuals$mosquito,
    individuals$human,
    states,
    variables,
    events$mosquito_infection
  )
  api <- mock_api(
    list(
      human = list(
        U = c(1, 2),
        A = c(3),
        D = c(4),
        age = c(20, 24, 5, 39),
        xi = c(.2, .3, .5, .9),
        ID = c(.2, .3, .5, .9)
      ),
      mosquito = list(
        Sm = c(1, 2, 3, 4),
        variety = c(1, 2, 3, 3)
      )
    ),
    parameters = parameters
  )

  mockery::stub(
    mosquito_infection_process,
    'bernoulli',
    mockery::mock(c(TRUE, TRUE, TRUE, FALSE))
  )

  mosquito_infection_process(api)

  updates <- mockery::mock_args(api$queue_state_update)

  expect_equal(updates[[1]][[1]]$name, 'mosquito')
  expect_equal(updates[[1]][[2]]$name, 'Pm')
  expect_equal(updates[[1]][[3]], c(1, 2, 3))
})
