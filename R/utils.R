vnapply <- function(X, FUN, ...) vapply(X, FUN, ..., numeric(1))

bernoulli <- function(size, p) runif(size, 0, 1) < p
