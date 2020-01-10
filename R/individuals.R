# Define population sizes
human_population <- 100 * 1000
mosquito_population <- 100 * human_population

create_states <- function() {
  list(
    # Human states
    S       = State$new("S", human_population),
    I       = State$new("I", 0),
    Treated = State$new("T", 0),
    D       = State$new("D", 0),
    A       = State$new("A", 0),
    U       = State$new("U", 0),
    # Mosquito states
    E       = State$new("E", mosquito_population),
    L       = State$new("L", 0),
    P       = State$new("P", 0),
    Sm      = State$new("Sm", 0),
    Im      = State$new("Im", 0)
  )
}

create_variables <- function(parameters) {
  initial_age <- trunc(rexp(human_population, rate=1/10))

  # Define variables
  age <- Variable$new("age", function(size) initial_age)
  last_bitten <- Variable$new("last_bitten", function(size) { rep(-1, size) })
  last_infected <- Variable$new("last_infected", function(size) { rep(-1, size) })
  infection_schedule <- Variable$new(
    "infection_schedule",
    function(size) { rep(-1, size) }
  )
  asymptomatic_infection_schedule <- Variable$new(
    "asymptomatic_infection_schedule",
    function(size) { rep(-1, size) }
  )

  # Maternal immunity
  icm <- Variable$new(
    "ICM",
    function(size) {
      first_immunity <- 1
      t <- initial_age * 365 * parameters$timestep_to_day
      first_immunity * exp(-(t * parameters$rm))
    }
  )

  # Pre-erythoctic immunity
  ib  <- Variable$new("IB", function(size) { rep(0, size) })
  # Acquired immunity to severe disease
  ica <- Variable$new("ICA", function(size) { rep(0, size) })

  xi <- Constant$new(
    "xi",
    function(n) {
      rlnorm(n, -parameters$sigma_squared/2,parameters$sigma_squared)
    }
  )

  mosquito_variety <- Constant$new(
    "variety",
    function(n) {
      p <- runif(n)
      v <- rep(0, n)
      v[which(p > .5)] <- 1
      v[which(p > .2 & p < .5)] <- 2
      v[which(p < .2)] <- 3
      v
    }
  )
  list(
    age = age,
    last_bitten = last_bitten,
    last_infected = last_infected,
    icm = icm,
    ib = ib,
    ica = ica,
    xi = xi,
    mosquito_variety = mosquito_variety,
    infection_schedule = infection_schedule,
    asymptomatic_infection_schedule = asymptomatic_infection_schedule
  )
}

create_individuals <- function(states, variables) {
  env <- environment()
  list2env(states, env)
  list2env(variables, env)

  human <- Individual$new(
    'human',
    list(S, I, Treated, D, A, U),
    variables = list(
      age,
      last_bitten,
      last_infected,
      ib,
      ica,
      icm,
      infection_schedule,
      asymptomatic_infection_schedule
    ),
    constants = list(xi)
  )

  mosquito <- Individual$new(
    'mosquito',
    list(E, L, P, Sm, Im),
    constants = list(mosquito_variety)
  )

  list(human = human, mosquito = mosquito)
}
