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
