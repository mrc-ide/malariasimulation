
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
        ib = c(.2, .3, .5, .9),
        xi = c(.2, .3, .5, .9),
        ica = c(.2, .3, .5, .9),
        iva = c(.2, .3, .5, .9),
        icm = c(.2, .3, .5, .9),
        ivm = c(.2, .3, .5, .9),
        infection_schedule = c(-1, 6, 3, 7),
        asymptomatic_infection_schedule = c(5, -1, -1, -1),
        last_infected = c(-1, -1, 1, -1),
        id = c(.2, .3, .5, .9)
      ),
      mosquito = list(
        Im = 1:100,
        mosquito_variety = c(rep(1, 25), rep(2, 25), rep(3, 50))
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
  updates <- infection_process(simulation_frame, 1, parameters)
  for (update in updates) {
    print(update$variable$name)
    print(update$index)
  }
})
