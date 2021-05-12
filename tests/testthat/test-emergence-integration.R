test_that('emergence process fails when there are not enough individuals', {
  parameters <- get_parameters()
  state <- individual::CategoricalVariable$new(
    c('Sm', 'Pm', 'Im', 'NonExistent'),
    c(rep('Im', 1000), rep('Sm', 1000))
  )
  species <- individual::CategoricalVariable$new(
    c('All'),
    rep('All', 2000)
  )
  emergence_process <- create_mosquito_emergence_process(
    list(),
    state,
    species,
    c('All'),
    parameters$dpl
  )
  mockery::stub(
    emergence_process,
    'solver_get_states', 
    mockery::mock(
      c(1000, 500, 100),
      c(1000, 500, 100)
    )
  )
  expect_error(emergence_process(0), '*')
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
  emergence_process <- create_mosquito_emergence_process(
    list(mockery::mock(), mockery::mock()),
    state,
    species,
    c('a', 'b'),
    parameters$dpl
  )

  mockery::stub(
    emergence_process,
    'solver_get_states', 
    mockery::mock(
      c(100000, 50000, 10000),
      c(10000, 5000, 1000)
    )
  )

  emergence_process(0)

  expect_bitset_update(
    state$queue_update,
    'Sm',
    seq(7777) + 2000
  )
  expect_bitset_update(
    species$queue_update,
    'a',
    seq(7777) + 2000
  )
  expect_bitset_update(
    state$queue_update,
    'Sm',
    seq(778) + 2000 + 7777,
    call = 2
  )
  expect_bitset_update(
    species$queue_update,
    'b',
    seq(778) + 2000 + 7777,
    call = 2
  )
})
