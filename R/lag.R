LaggedValue <- R6::R6Class(
  'LaggedValue',
  private = list(
    history = NULL
  ),
  public = list(
    initialize = function(max_lag, default) {
      private$history <- create_timeseries(max_lag, default)
    },

    save = function(value, timestep) {
      timeseries_push(private$history, value, timestep)
    },
    
    get = function(timestep) {
      timeseries_at(private$history, timestep, TRUE)
    },

    save_state = function() {
      timeseries_save_state(private$history)
    },

    restore_state = function(t, state) {
      timeseries_restore_state(private$history, state)
    }
  )
)
