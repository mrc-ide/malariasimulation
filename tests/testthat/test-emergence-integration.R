
test_that('egg_laying_process fails when there are not enough individuals', {
  skip("to be translated to cpp")
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events)
  egg_laying_process <- create_egg_laying_process(
    individuals$mosquito,
    states$Sm,
    states$Im,
    states$Unborn,
    states$E
  )

  api <- mock_api(
    list(
      mosquito = list(
        Sm = seq_len(1000),
        Im = seq_len(1000) + 1000,
        Unborn = c()
      )
    )
  )
  expect_error(
    egg_laying_process(api),
    '*'
  )

  api <- mock_api(
    list(
      mosquito = list(
        Sm = seq_len(1000),
        Im = seq_len(1000) + 1000,
        Unborn = seq_len(1000) + 2000
      )
    )
  )
  expect_error(
    egg_laying_process(api),
    '*'
  )
})

test_that('egg_laying_process creates the correct number of larvae', {
  skip("to be translated to cpp")
  parameters <- get_parameters()
  parameters$beta <- 5
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events)
  egg_laying_process <- create_egg_laying_process(
    individuals$mosquito,
    states$Sm,
    states$Im,
    states$Unborn,
    states$E
  )

  api <- mock_api(
    list(
      mosquito = list(
        Sm = seq_len(1000),
        Im = seq_len(1000) + 1000,
        Unborn = seq_len(100000) + 2000
      )
    ),
    parameters = parameters
  )
  egg_laying_process(api)
  expect_equal(update$individual$name, 'mosquito')
  expect_equal(update$state$name, 'E')
  expect_equal(length(update$index), 10000)
  expect(all(update$index >= 2000 & update$index < 100000 + 2000), 'incorrect range')
})

test_that('larval_death_process works with no larvae', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events)
  larval_death_process <- create_larval_death_process(
    individuals$mosquito,
    states$E,
    states$L,
    states$Unborn,
    events
  )
  api <- mock_api(
    list(
      mosquito = list(
        E = c(),
        L = c(),
        Unborn = seq_len(100)
      )
    ),
    timestep = 100,
    parameters = parameters
  )
  larval_death_process(api)
  updates <- mockery::mock_args(api$queue_state_update)
  expect_equal(updates[[1]][[1]]$name, 'mosquito')
  expect_equal(updates[[1]][[2]]$name, 'Unborn')
  expect_length(updates[[1]][[3]], 0)
})

test_that('larval_death_process kills the expected larvae', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events)
  larval_death_process <- create_larval_death_process(
    individuals$mosquito,
    states$E,
    states$L,
    states$Unborn,
    events
  )
  api <- mock_api(
    list(
      mosquito = list(
        E = c(1, 2, 3, 4),
        L = c(5, 6, 7, 8),
        Unborn = seq_len(100) + 8
      )
    ),
    timestep = 100,
    parameters = parameters
  )

  mocked <- mockery::stub(
    larval_death_process,
    'bernoulli',
    mockery::mock(
      c(TRUE, TRUE, FALSE, FALSE),
      c(TRUE, FALSE, TRUE, FALSE)
    )
  )

  larval_death_process(api)
  updates <- mockery::mock_args(api$queue_state_update)
  expect_equal(updates[[1]][[1]]$name, 'mosquito')
  expect_equal(updates[[1]][[2]]$name, 'Unborn')
  expect_equal(updates[[1]][[3]], c(1, 2, 5, 7))
})
