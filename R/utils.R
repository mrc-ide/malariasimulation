vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, ..., numeric(1))
}
