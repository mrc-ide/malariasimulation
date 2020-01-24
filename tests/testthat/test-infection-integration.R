
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

  epsilon <- c(0.2872749, 0.6559928, 1.6966581, 6.2333683)
  b <- c(0.5900733, 0.59006960, 0.59006960)
  phi <- c(0.07497956, 0.07496497)
  theta <- c(0.07360586, 0.07194696)
  mockery::stub(
    infection_process,
    'runif',
    mock_returns(list(
      mock_random(epsilon, c(TRUE, TRUE, TRUE, FALSE)),
      mock_random(b, c(TRUE, TRUE, FALSE)),
      mock_random(phi, c(TRUE, FALSE)),
      mock_random(theta, c(FALSE, TRUE))
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
