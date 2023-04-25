#' Create progress process
#'
#' @param timesteps Simulation timesteps
#'
#' @return Progress bar process function
create_progress_process <- function(timesteps){
  pb <- progress::progress_bar$new(
    format = "  running [:bar] :percent eta: :eta",
    total = timesteps, clear = FALSE, width= 60)
  
  function(timestep) {
    pb$tick()
  }
}