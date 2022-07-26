
test_that('Asymptomatic infection listener updates infectivity correctly', {
  parameters <- get_parameters(list(human_population = 4))
  variables <- create_variables(parameters)

  update_mock <- mockery::mock()
  variables$infectivity <- list(queue_update = update_mock)

  to_move <- individual::Bitset$new(4)$insert(c(2, 4))
  mockery::stub(
    update_to_asymptomatic_infection,
    'asymptomatic_infectivity',
    mockery::mock(c(.2, .5))
  )
  update_to_asymptomatic_infection(variables, parameters, 1, to_move)

  expect_bitset_update(
    update_mock,
    c(.2, .5),
    c(2, 4)
  )
})
