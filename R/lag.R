LaggedValue <- R6::R6Class(
  'LaggedValue',
  private = list(
    history = NULL
  ),
  public = list(
    initialize = function(max_lag, default) {
      private$history <- create_history(max_lag, default)
    },

    save = function(value, timestep) {
      history_push(private$history, value, timestep)
    },
    
    get = function(timestep) {
      history_at(private$history, timestep)
    }
  )
)
