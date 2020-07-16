vnapply <- function(X, FUN, ...) vapply(X, FUN, ..., numeric(1))

bernoulli <- function(size, p) runif(size, 0, 1) < p

discretise_normal <- function(values, n_groups) {
  quads <- statmod::gauss.quad.prob(n_groups, dist='normal')
  as.numeric(cut(
    values,
    quads$nodes,
    include.lowest = TRUE,
    right = FALSE
  ))
}

get_age <- function(birth_timesteps, current_timestep) {
  age <- current_timestep - birth_timesteps
  age[age >= 100 * 365] <- 99 *365
  age
}

remove_keys <- function(x, n) { for (name in n) { x[[name]] <- NULL }; x }
