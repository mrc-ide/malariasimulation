#' @title Mosquito births
#' @description
#' This is the process for mosquito birth, it defines how many new early stage
#' larvae are created on each timestep.
#' @param mosquito, the mosquito individual
#' @param Sm, the susceptable mosquito state
#' @param Im, the infected mosquito state
#' @param Unborn, the unborn mosquito state
#' @param E, the early stage larval state
create_egg_laying_process <- function(mosquito, Sm, Im, Unborn, E) {
  function(simulation_frame, timestep, parameters) {
    m <- simulation_frame$get_state(mosquito, Sm, Im)
    unborn <- simulation_frame$get_state(mosquito, Unborn)
    if (length(m) > 0) {
      n_eggs <- parameters$beta * length(m)
      if (n_eggs > length(unborn)) {
        stop('Run out of mosquitos')
      }
      if (n_eggs >= 1) {
        return(individual::StateUpdate$new(mosquito, E, unborn[seq_len(n_eggs)]))
      }
    }
  }
}
