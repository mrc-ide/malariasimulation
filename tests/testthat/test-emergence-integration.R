test_that('emergence process fails when there are not enough individuals', {
  parameters <- get_parameters()
  state <- individual::CategoricalVariable$new(
    c('Sm', 'Pm', 'Im', 'NonExistent'),
    c(rep('Im', 1000), rep('Sm', 1000))
  )
  species <- individual::CategoricalVariable$new(
    c('gamb'),
    rep('gamb', 2000)
  )
  solvers <- list(
    mock_solver(c(1000, 500, 100)),
    mock_solver(c(1000, 500, 100))
  )

  emergence_process <- create_mosquito_emergence_process(
    solvers,
    state,
    species,
    c('gamb'),
    parameters$dpl
  )
  expect_error(emergence_process(0), 'Not enough mosquitoes')
})

test_that('emergence_process creates the correct number of susceptables', {
  parameters <- get_parameters()
  state <- mock_category(
    c('Sm', 'Pm', 'Im', 'NonExistent'),
    c(rep('Im', 1000), rep('Sm', 1000), rep('NonExistent', 100000))
  )
  species <- mock_category(
    c('a', 'b'),
    c('a', 'b')
  )
  solvers <- list(
    mock_solver(c(100000, 50000, 10000)),
    mock_solver(c(1000, 5000, 1000))
  )

  emergence_process <- create_mosquito_emergence_process(
    solvers,
    state,
    species,
    c('a', 'b'),
    parameters$dpl
  )

  emergence_process(0)

  expect_bitset_update(
    state$queue_update_mock(),
    'Sm',
    seq(7777) + 2000
  )
  expect_bitset_update(
    species$queue_update_mock(),
    'a',
    seq(7777) + 2000
  )
  expect_bitset_update(
    state$queue_update_mock(),
    'Sm',
    seq(778) + 2000 + 7777,
    call = 2
  )
  expect_bitset_update(
    species$queue_update_mock(),
    'b',
    seq(778) + 2000 + 7777,
    call = 2
  )
})
