
test_that('Asymptomatic infection listener updates infectivity correctly', {
  parameters <- get_parameters(list(human_population = 4))
  variables <- create_variables(parameters)
  listener <- create_asymptomatic_update_listener(variables, parameters)

  mockery::stub(listener, 'asymptomatic_infectivity', mockery::mock(c(.2, .5)))

  update_mock <- mockery::mock()
  variables$infectivity <- list(queue_update = update_mock)

  to_move <- individual::Bitset$new(4)
  to_move$insert(c(2, 4))
  listener(1, to_move)

  expect_bitset_update(
    update_mock,
    c(.2, .5),
    c(2, 4)
  )
})
