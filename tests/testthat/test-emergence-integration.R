
test_that('egg_laying_process fails when there are not enough individuals', {
  parameters <- get_parameters()
  bind_process_to_default_model(egg_laying_process, parameters)
  simulation_frame <- mock_simulation_frame(
    list(
      mosquito = list(
        Sm = seq_len(1000),
        Im = seq_len(1000) + 1000,
        Unborn = c()
      )
    )
  )
  expect_error(
    egg_laying_process(simulation_frame, 1, parameters),
    '*'
  )

  simulation_frame <- mock_simulation_frame(
    list(
      mosquito = list(
        Sm = seq_len(1000),
        Im = seq_len(1000) + 1000,
        Unborn = seq_len(1000) + 2000
      )
    )
  )
  expect_error(
    egg_laying_process(simulation_frame, 1, parameters),
    '*'
  )
})

test_that('egg_laying_process creates the correct number of larvae', {
  parameters <- get_parameters()
  parameters$beta <- 5
  bind_process_to_default_model(egg_laying_process, parameters)
  simulation_frame <- mock_simulation_frame(
    list(
      mosquito = list(
        Sm = seq_len(1000),
        Im = seq_len(1000) + 1000,
        Unborn = seq_len(100000) + 2000
      )
    )
  )
  update <- egg_laying_process(simulation_frame, 1, parameters)
  expect_equal(update$individual$name, 'mosquito')
  expect_equal(update$state$name, 'E')
  expect_equal(length(update$index), 10000)
  expect(all(update$index >= 2000 && update$index < 100000 + 2000), 'incorrect range')
})

test_that('larval_death_process works with no larvae', {
})

test_that('larval_death_process kills the expected larvae', {
})
