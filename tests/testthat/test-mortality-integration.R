test_that('mortality_process resets humans correctly', {
  parameters <- get_parameters()
  bind_process_to_default_model(mortality_process, parameters)
  simulation_frame <- mock_simulation_frame(
    list(
      human = list(
        Treated = c(1),
        D = c(2),
        age = c(20, 24, 5, 39),
        is_severe = c(1., 1., 0., 0.),
        xi_group = c(1, 1, 2, 2),
        ICM = c(1, 2, 3, 4),
        IVM = c(1, 2, 3, 4)
      )
    )
  )

  mockery::stub(
    mortality_process,
    'bernoulli',
    mock_returns(list(
      c(FALSE, FALSE, FALSE, TRUE),
      c(TRUE),
      c(FALSE, TRUE)
    )),
    depth = 2
  )

  mockery::stub(
    mortality_process,
    'sample',
    mock_returns(list(c(1), c(4))),
    depth = 2
  )

  updates <- mortality_process(simulation_frame, 1, parameters)

  died <- c(2, 4)

  expect_any(updates, function(update) {
    all(
      update$variable$name == 'age',
      update$value == 0,
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'last_bitten',
      update$value == -1,
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'last_infected',
      update$value == -1,
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'infection_schedule',
      update$value == -1,
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'asymptomatic_infection_schedule',
      update$value == -1,
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'ICM',
      update$value == c(2 * parameters$pm, 4 * parameters$pm),
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'IVM',
      update$value == c(2 * parameters$pm, 4 * parameters$pm),
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'IB',
      update$value == -1,
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'ICA',
      update$value == -1,
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'IVA',
      update$value == -1,
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'ID',
      update$value == -1,
      setequal(update$index, died)
    )
  })
  expect_any(updates, function(update) {
    all(
      update$variable$name == 'is_severe',
      update$value == 0,
      setequal(update$index, died)
    )
  })
})
