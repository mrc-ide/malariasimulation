vnapply <- function(X, FUN, ...) vapply(X, FUN, ..., numeric(1))

uniform_gt <- function(size, p) runif(size, 0, 1) > p
uniform_lt <- function(size, p) runif(size, 0, 1) < p
