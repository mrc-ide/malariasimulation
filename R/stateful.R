#' @title Save the state of a list of \emph{stateful objects}
#' @description The state of each element is saved and stored into a single
#' object, representing them in a way that can be exported and re-used later.
#' @param objects A list of stateful objects to be saved. Stateful objects are
#' instances of R6 classes with a pair of \code{save_state} and
#' \code{restore_state} methods.
#' @noRd
save_state <- function(objects) {
  lapply(objects, function(o) o$save_state())
}

#' @title Restore the state of a collection of stateful objects
#' @description This is the counterpart of \code{save_state}. Calling it
#' restores the collection of objects back into their original state.
#' @param state A state object, as returned by the \code{save_state} function.
#' @param objects A collection of stateful objects to be restored.
#' @noRd
restore_state <- function(state, objects) {
  stopifnot(length(state) == length(objects))
  for (i in seq_along(state)) {
    objects[[i]]$restore_state(state[[i]])
  }
}

#' @title a placeholder class to save the random number generator class.
#' @description the class integrates with the simulation loop to save and
#' restore the random number generator class when appropriate.
#' @noRd
RandomState <- R6::R6Class(
  'RandomState',
  public = list(
    save_state = function() {
      random_save_state()
    },
    restore_state = function(state) {
      random_restore_state(state)
    }
  )
)
