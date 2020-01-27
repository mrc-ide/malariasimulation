
test_that('human infection_process creates the correct updates', {
  parameters <- get_parameters()
  bind_process_to_default_model(infection_process, parameters)
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
    'uniform_gt',
    mock_returns(list(
      c(TRUE, TRUE, TRUE, FALSE),
      c(TRUE, TRUE, FALSE),
      c(TRUE, FALSE),
      c(FALSE, TRUE)
    ))
  )
  updates <- infection_process(simulation_frame, 5, parameters)
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'ICA',
      all.equal(update$value, c(1.2, 1.3)) == TRUE,
      all.equal(update$index, c(1, 2)) == TRUE
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'IVA',
      all.equal(update$value, c(1.2, 1.3)) == TRUE,
      all.equal(update$index, c(1, 2)) == TRUE
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'ID',
      all.equal(update$value, c(1.2, 1.3)) == TRUE,
      all.equal(update$index, c(1, 2)) == TRUE
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'IB',
      all.equal(update$value, c(1.2, .3, .5)) == TRUE,
      all.equal(update$index, c(1, 2, 3)) == TRUE
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'infection_schedule',
      update$value == 5 + parameters$de,
      all.equal(update$index, c(1, 2)) == TRUE
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'asymptomatic_infection_schedule',
      all.equal(update$index, c()) == TRUE
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'last_bitten',
      update$value == 5,
      all.equal(update$index, c(1, 2, 3))
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'last_infected',
      update$value == 5,
      all.equal(update$index, c(1, 2))
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'is_severe',
      update$value == 1,
      all.equal(update$index, c(2))
    )
  })
})

test_that('mosquito_infection_process creates the correct updates', {
  parameters <- get_parameters()
  bind_process_to_default_model(mosquito_infection_process, parameters)
  simulation_frame <- mock_simulation_frame(
    list(
      human = list(
        Treated = c(1),
        U = c(2),
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
    'uniform_gt',
    mock_returns(list(
      c(TRUE, TRUE, TRUE, FALSE)
    ))
  )

  update <- mosquito_infection_process(simulation_frame, 5, parameters)

  expect_equal(update$individual$name, 'mosquito')
  expect_equal(update$state$name, 'Im')
  expect_equal(update$index, c(1, 2, 3))
})
