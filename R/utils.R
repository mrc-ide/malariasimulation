vnapply <- function(X, FUN, ...) vapply(X, FUN, ..., numeric(1))

bernoulli <- function(size, p) runif(size, 0, 1) < p

discretise <- function(values, n_groups) {
  as.numeric(cut(
    values,
    seq(min(values), max(values), length.out=n_groups + 1),
    include.lowest = TRUE,
    right = FALSE,
    labels = seq_len(n_groups)
  ))
}
