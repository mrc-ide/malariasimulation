#' @title GridRender
#' @description Class to render grid based output for the simulation.
#' @importFrom R6 R6Class
#' @export
GridRender <- R6::R6Class(
  'GridRender',
  private = list(
    vectors = list(),
    timesteps = 0,
    width = 0
  ),
  public = list(

    #' @description
    #' Initialise a grid renderer for the simulation, creates the default state
    #' renderers.
    #' @param timesteps number of timesteps in the simulation.
    #' @param width of the grid
    initialize = function(timesteps, width) {
      private$timesteps = timesteps
      private$width = width
      private$vectors[['timestep']] <- seq_len(timesteps)
    },
    
    #' @description
    #' Set a default value for a rendered output
    #' renderers.
    #' @param name the variable to set a default for.
    #' @param value  the default value to set for a variable.
    set_default = function(name, value) {
      if (name == 'timestep') {
        stop("Cannot set default value for variable 'timestep'")
      }
      private$vectors[[name]] = matrix(
        value,
        nrow = private$timesteps,
        ncol = private$width
      )
    },

    #' @description
    #' Update the render with new simulation data.
    #' @param name the variable to render.
    #' @param value the value to store for the variable.
    #' @param timestep the time-step of the data point.
    render = function(name, value, timestep) {
      if (name == 'timestep') {
        stop("Please don't name your variable 'timestep'")
      }
      if (!(name %in% names(private$vectors))) {
        private$vectors[[name]] = matrix(
          NA,
          nrow = private$timesteps,
          ncol = private$width
        )
      }
      private$vectors[[name]][timestep,] = value
    },

    #' @description
    #' Return the render as a list of matrices
    get_values = function() {
      private$vectors
    }
  )
)

#' @title Render incidence statistics
#' 
#' @description renders incidence (new for this timestep) for indivduals
#' 
#' @param birth variable for birth of the individual
#' @param renderer object for model outputs
#' @param target incidence population
#' @param source_pop the population which is sampled for infection
#' @param prob probability of infection
#' @param prefix for model outputs
#' @param lowers age bounds
#' @param uppers age bounds
#' @param timestep current target
#' 
#' @noRd
incidence_grid_renderer <- function(
  birth,
  renderer,
  target,
  prefix,
  timestep
  ) {
  renderer$render(
    paste0('n_', prefix),
    grid_count(birth, target, timestep),
    timestep
  )
}

grid_count <- function(birth, selected, timestep) {
  if (is.null(selected)) {
    selected_births <- birth$get_values()
  } else {
    selected_births <- birth$get_values(selected)
  }
  age <- floor(get_age(selected_births, timestep) / 365)
  age[age < 0] <- NA
  age[age > 100] <- NA 
  non_zero <- table(age)
  counts <- rep(0, 101)
  counts[as.numeric(names(non_zero)) + 1] <- non_zero
  counts
}
