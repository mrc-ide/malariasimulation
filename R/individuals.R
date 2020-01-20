# Define population sizes
human_population <- 100 * 1000
mosquito_limit <- 100 * human_population
n_heterogeneity_groups <- 5

create_states <- function() {
  list(
    # Human states
    S       = individual::State$new("S", human_population),
    I       = individual::State$new("I", 0),
    Treated = individual::State$new("T", 0),
    D       = individual::State$new("D", 0),
    A       = individual::State$new("A", 0),
    U       = individual::State$new("U", 0),
    # Mosquito states
    E       = individual::State$new("E", 0),
    L       = individual::State$new("L", 0),
    P       = individual::State$new("P", 0),
    Sm      = individual::State$new("Sm", mosquito_limit %/% 2),
    Im      = individual::State$new("Im", 0),
    Unborn  = individual::State$new("Unborn", mosquito_limit %/% 2)
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

  ivm <- Variable$new(
    "IVM",
    function(size) {
      first_immunity <- 1
      t <- initial_age * 365 * parameters$timestep_to_day
      first_immunity * exp(-(t * parameters$rm))
    }
  )

  # Pre-erythoctic immunity
  ib  <- Variable$new("IB", function(size) { rep(0, size) })
  # Acquired immunity to clinical disease
  ica <- Variable$new("ICA", function(size) { rep(0, size) })
  # Acquired immunity to severe disease
  iva <- Variable$new("IVA", function(size) { rep(0, size) })
  # Acquired immunity to detectability
  id <- Variable$new("ID", function(size) { rep(0, size) })

  # is_severe, 1 iff the individual is currently affected by severe disease
  is_severe <- Variable$new("is_severe", function(size) { rep(0, size) })

  xi_values <- rlnorm(n, -parameters$sigma_squared/2,parameters$sigma_squared)
  xi <- Constant$new(
    "xi",
    function(n) {
      xi_values
    }
  )

  xi_group <- Constant$new(
    "xi_group",
    function(n) {
      discretise(xi_values, n_heterogeneity_groups)
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
    iva = iva,
    id = id,
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
      iva,
      id,
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

discretise <- function(values, n_groups) {
  as.numeric(cut(
    values,
    seq(min(values), max(values), length.out=n_groups + 1),
    include.lowest = TRUE,
    right = FALSE,
    labels = seq_len(n_groups)
  ))
}
