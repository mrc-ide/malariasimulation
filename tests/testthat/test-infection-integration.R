
test_that('human infection_process creates the correct updates', {
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables)

  infection_process <- create_infection_process(
    individuals$human,
    individuals$mosquito,
    states,
    variables
  )

  simulation_frame <- mock_simulation_frame(
    list(
      human = list(
        S = c(1),
        U = c(2),
        A = c(3),
        D = c(4),
        age = c(20, 24, 5, 39),
        IB = c(.2, .3, .5, .9),
        xi = c(.2, .3, .5, .9),
        ICA = c(.2, .3, .5, .9),
        IVA = c(.2, .3, .5, .9),
        ICM = c(.2, .3, .5, .9),
        IVM = c(.2, .3, .5, .9),
        infection_schedule = c(-1, 1, 3, 7),
        asymptomatic_infection_schedule = c(5, -1, -1, -1),
        last_infected = c(-1, -1, 1, -1),
        last_bitten = c(-1, 1, 1, -1),
        ID = c(.2, .3, .5, .9)
      ),
      mosquito = list(
        Im = 1:100,
        variety = c(rep(1, 25), rep(2, 25), rep(3, 50))
      )
    )
  )

  mockery::stub(
    infection_process,
    'bernoulli',
    mockery::mock(
      c(FALSE, FALSE, TRUE),
      c(TRUE, FALSE),
      c(FALSE, TRUE)
    )
  )

  updates <- infection_process(simulation_frame, 5, parameters)

  expect_any(updates, function(update) {
    all(
      update$variable$name == 'ICA',
      setequal(update$value, c(1.3, 0.5)),
      setequal(update$index, c(2, 3))
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'IVA',
      setequal(update$value, c(1.3, 0.5)),
      setequal(update$index, c(2, 3))
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'ID',
      setequal(update$value, c(1.3, 0.5)),
      setequal(update$index, c(2, 3))
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'IB',
      setequal(update$value, c(.3, .5, 1.9)),
      setequal(update$index, c(2, 3, 4))
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'infection_schedule',
      update$value == 5 + parameters$de,
      setequal(update$index, c(2, 3))
    )
  })
  expect_none(updates, function(update) {
    update$variable$name == 'asymptomatic_infection_schedule'
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'last_bitten',
      update$value == 5,
      setequal(update$index, c(2, 3, 4))
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'last_infected',
      update$value == 5,
      setequal(update$index, c(2, 3))
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'is_severe',
      update$value == 1,
      setequal(update$index, c(3))
    )
  })
})

test_that('mosquito_infection_process creates the correct updates', {
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables)
  mosquito_infection_process <- create_mosquito_infection_process(
    individuals$mosquito,
    individuals$human,
    states,
    variables
  )
  simulation_frame <- mock_simulation_frame(
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
    )
  )

  mockery::stub(
    mosquito_infection_process,
    'bernoulli',
    mockery::mock(c(TRUE, TRUE, TRUE, FALSE))
  )

  update <- mosquito_infection_process(simulation_frame, 5, parameters)

  expect_equal(update$individual$name, 'mosquito')
  expect_equal(update$state$name, 'Im')
  expect_equal(update$index, c(1, 2, 3))
})
