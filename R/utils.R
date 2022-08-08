vnapply <- function(X, FUN, ...) vapply(X, FUN, ..., numeric(1))
vlapply <- function(X, FUN, ...) vapply(X, FUN, ..., logical(1))
vcapply <- function(X, FUN, ...) vapply(X, FUN, ..., character(1))

#' @importFrom stats rbinom 
bernoulli <- function(size, p) sample.int(size, rbinom(1, size, min(p, 1)))

sample_bitset <- function(b, rate) {
  individual::filter_bitset(b, bernoulli(b$size(), rate))
}

bitset_at <- function(b, i) {
  individual::filter_bitset(b, i)
}

bernoulli_multi_p <- function(p) bernoulli_multi_p_cpp(p)


#' @title find the indices of a where it intersects with b
#' @description synonymous with \code{which(a$to_vector() %in%
#' b$to_vector())} but faster
#' @param a the bitset to index
#' @param b the bitset to check
#' @noRd
bitset_index <- function(a, b) bitset_index_cpp(a$.bitset, b$.bitset)

#' @importFrom stats runif
log_uniform <- function(size, rate) -rate * log(runif(size))

approx_sum <- function(X, n) abs(sum(X) - n) < sqrt(.Machine$double.eps)

get_age <- function(birth_timesteps, current_timestep) {
  current_timestep - birth_timesteps
}

remove_keys <- function(x, n) { for (name in n) { x[[name]] <- NULL }; x }

invlogit <- function(x) exp(x) / (1 + exp(x))

to_char_vector <- function(v) vapply(v, function(n) toString(n), character(1))

weibull_survival <- function(t, shape, scale) exp(-((t/scale)^shape))

between <- function(x, l, u) (x >= l) && (x <= u)

#'@title Inverse truncated exponential CDF
#'@param u quantile
#'@param m rate
#'@param t truncation
#'@noRd
itexp <- function(u, m, t) -log(1 - u * (1 - exp(-t * m)))/m

#'@title Truncated exponential sample
#'@param n number of samples
#'@param m rate
#'@param t truncation
#'@noRd
rtexp <- function(n, m, t) { itexp(runif(n), m, t) }

#'@title Find index of the latest timestep in vector of timesteps
#'@param ts timesteps, assumed to be sorted and unique
#'@param t current timestep
#'@noRd
match_timestep <- function(ts, t) {
  min(sum(ts <= t), length(ts))
}
