## Define classes to resolve competing hazards
CompetingOutcome <- R6::R6Class(
  "CompetingOutcome",
  private = list(
    targeted_process = NULL
  ),
  public = list(
    initialize = function(targeted_process, size){
      if (!is.function(targeted_process)){
        stop("targeted_process must be a function")
      }
      if (!is.numeric(size) || size <= 0){
        stop("size must be positive integer")
      }
      private$targeted_process <- targeted_process
      self$rates <- rep(0, size)
    },
    set_rates = function(rates){
      self$rates <- rates
    },
    execute = function(t, target){
      private$targeted_process(t, target)
      self$rates <- rep(0, length(self$rates))
    },
    rates = NULL
  )
)

CompetingHazard <- R6::R6Class(
  "CompetingHazard",
  private = list(
    outcomes = list(),
    size = NULL,
    # RNG is passed in because mockery is not able to stub runif
    # TODO: change when fixed
    rng = NULL
  ),
  public = list(
    initialize = function(outcomes, rng = runif){
      if (length(outcomes) == 0){
        stop("At least one outcome must be provided")
      }
      if (!all(sapply(outcomes, function(x) inherits(x, "CompetingOutcome")))){
        stop("All outcomes must be of class CompetingOutcome")
      }
      private$outcomes <- outcomes
      private$size <- length(outcomes[[1]]$rates)
      private$rng <- rng
    },
    resolve = function(t){
      event_rates <- do.call(
        'cbind',
        lapply(private$outcomes, function(x) x$rates)
      )
      occur_rates <- rowSums(event_rates)
      occur_rng <- private$rng(private$size)
      occurs <- occur_rng < rate_to_prob(occur_rates)
      norm_probs <- event_rates / occur_rates
      norm_probs[is.na(norm_probs)] <- 0
      cum_probs <- apply(norm_probs, 1, cumsum)
      event_rng <- private$rng(private$size)
      for(o in seq_along(private$outcomes)){
        if (o == 1) {
          selected <- (event_rng <= cum_probs[o,])
        } else {
          selected <- (event_rng > cum_probs[o - 1,]) &
            (event_rng <= cum_probs[o,])
        }
        target <- individual::Bitset$new(private$size)$insert(
          which(selected & occurs)
        )
        if (target$size() > 0){
          private$outcomes[[o]]$execute(t, target)
        }
      }
    }
  )
)
