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

      self$target <- individual::Bitset$new(size)
      self$rates <- NULL
    },
    set_rates = function(target, rates){
      stopifnot(target$size() == length(rates))

      self$target$copy_from(target)
      self$rates <- rates
    },
    execute = function(t, target){
      private$targeted_process(t, target)
    },
    reset = function() {
      self$target$clear()
      self$rates <- NULL
    },
    target = NULL,
    rates = NULL
  )
)

CompetingHazard <- R6::R6Class(
  "CompetingHazard",
  private = list(
    size = NULL,
    outcomes = list(),
    # RNG is passed in because mockery is not able to stub runif
    # TODO: change when fixed
    rng = NULL
  ),
  public = list(
    initialize = function(size, outcomes, rng = runif){
      if (length(outcomes) == 0){
        stop("At least one outcome must be provided")
      }
      if (!all(sapply(outcomes, function(x) inherits(x, "CompetingOutcome")))){
        stop("All outcomes must be of class CompetingOutcome")
      }
      private$size <- size
      private$outcomes <- outcomes
      private$rng <- rng
    },
    resolve = function(t){
      candidates <- individual::Bitset$new(private$size)
      for (o in private$outcomes) {
        candidates$or(o$target)
      }

      rates <- matrix(ncol = length(private$outcomes), nrow = candidates$size(), 0)
      for (i in seq_along(private$outcomes)) {
        idx <- bitset_index(
          candidates,
          private$outcomes[[i]]$target)

        rates[idx, i] <- private$outcomes[[i]]$rates
      }

      total_rates <- rowSums(rates)
      probs <- rate_to_prob(total_rates) * (rates / total_rates)
      probs[is.na(probs)] <- 0

      rng <- private$rng(candidates$size())

      cumulative <- rep(0, candidates$size())

      for (o in seq_along(private$outcomes)) {
        next_cumulative <- cumulative + probs[,o]
        selected <- (rng > cumulative) & (rng <= next_cumulative)
        cumulative <- next_cumulative

        target <- bitset_at(candidates, selected)
        if (target$size() > 0) {
          private$outcomes[[o]]$execute(t, target)
        }
        private$outcomes[[o]]$reset()
      }
    }
  )
)
