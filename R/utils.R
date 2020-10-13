vnapply <- function(X, FUN, ...) vapply(X, FUN, ..., numeric(1))
vlapply <- function(X, FUN, ...) vapply(X, FUN, ..., logical(1))

#' @importFrom stats rbinom 
bernoulli <- function(size, p) sample.int(size, rbinom(1, size, min(p, 1)))

#' @importFrom stats runif
bernoulli_multi_p <- function(size, p) runif(size, 0, 1) < p

#' @importFrom stats runif
log_uniform <- function(size, rate) -rate * log(runif(size))

approx_sum <- function(X, n) abs(sum(X) - n) < sqrt(.Machine$double.eps)

discretise_normal <- function(values, n_groups) {
  quads <- statmod::gauss.quad.prob(n_groups, dist='normal')
  breaks <- quads$nodes
  breaks[[1]] <- min(min(values), min(breaks))
  breaks <- c(breaks, max(max(values), max(breaks) + 1))
  as.numeric(cut(
    values,
    breaks,
    include.lowest = TRUE,
    right = FALSE
  ))
}

get_age <- function(birth_timesteps, current_timestep) {
  current_timestep - birth_timesteps
}

remove_keys <- function(x, n) { for (name in n) { x[[name]] <- NULL }; x }

invlogit <- function(x) exp(x) / (1 + exp(x))

small_delta <- 1e-16
