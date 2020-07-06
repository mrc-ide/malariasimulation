create_prevelance_renderer <- function(human, D, A, birth) {
  function(api) {
    age <- get_age(api$get_variable(human, birth), api$get_timestep())
    in_range <- which((age > 2 * 365) & (age < 11 * 365))
    prev <- length(
      intersect(api$get_state(human, D), in_range)
    ) + length(
      intersect(api$get_state(human, A), in_range)
    )
    api$render('prev_2_10', prev)
    api$render('n_2_10', length(in_range))
  }
}

create_age_dist_renderer <- function(human, birth) {
  function(api) {
    age <- trunc(get_age(api$get_variable(human, birth), api$get_timestep()) / 365)
    api$render(
      'age_dist',
      vapply(seq(100) - 1, function(a) sum(a == age), numeric(1))
    )
  }
}

create_ica_dist_renderer <- function(human, ica, birth) {
  function(api) {
    age <- trunc(get_age(api$get_variable(human, birth), api$get_timestep()) / 365)
    ica_v <- api$get_variable(human, ica)
    api$render(
      'ica_dist',
      vapply(seq(100) - 1, function(a) mean(ica_v[a == age]), numeric(1))
    )
  }
}
