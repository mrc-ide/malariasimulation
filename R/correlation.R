INTS <- c(
  'rtss',
  'mda',
  'smc',
  'tbv',
  'bednets',
  'spraying'
)

#' Class: Correlation parameters
#' Describes an event in the simulation
#' @export CorrelationParameters
CorrelationParameters <- R6::R6Class(
  'CorrelationParameters',
  private = list(
    interventions = list(),
    n_ints = function() length(private$interventions),
    population = 0,
    rho_matrix = NULL,
    rho = function() diag(private$rho_matrix),
    .sigma = NULL,
    .mvnorm = NULL
  ),
  public = list(

    #' @description initialise correlation parameters
    #' @param parameters model parameters
    initialize = function(parameters) {
      # Find a list of enabled interventions
      enabled <- vlapply(INTS, function(name) parameters[[name]])
      private$interventions <- INTS[enabled]

      # Initialise a rho matrix for our interventions
      n_ints <- private$n_ints()
      private$rho_matrix <- matrix(
        0,
        nrow = n_ints,
        ncol = n_ints,
        dimnames = list(private$interventions, private$interventions)
      )

      # Store population for mvnorm draws
      private$population <- parameters$human_population
    },

    #' @description Add rho between rounds
    #' @param intervention string representing the intervention to update
    #' @param rho value between 0 and 1 representing the correlation between rounds of
    #' the intervention
    inter_round_rho = function(int, rho) {
      if (!(int %in% private$interventions)) {
        stop(paste0('invalid intervention name: ', int))
      }
      if (rho < 0 || rho > 1) {
        stop(paste0('rho for ', int, 'must be between 0 and 1'))
      }
      if (rho == 1) {
        rho <- 1 - .Machine$double.eps
      }
      private$rho_matrix[[int, int]] <- rho
      private$.sigma <- NULL
    },

    #' @description Add rho between interventions
    #' @param int_1 string representing the first intervention
    #' @param int_2 string representing the second intervention (intechangable
    #' with int_1)
    #' @param rho value between -1 and 1 representing the correlation between rounds of
    #' the intervention
    inter_intervention_rho = function(int_1, int_2, rho) {
      if (!(int_1 %in% private$interventions)) {
        stop(paste0('invalid intervention name: ', int))
      }
      if (!(int_2 %in% private$interventions)) {
        stop(paste0('invalid intervention name: ', int))
      }
      if (rho < -1 || rho > 1) {
        stop(paste0(
          'rho between ',
          int_1,
          ' and ',
          int_2,
          'must be between -1 and 1'
        ))
      }
      private$rho_matrix[[int_1, int_2]] <- rho
      private$rho_matrix[[int_2, int_1]] <- rho
      private$.sigma <- NULL
    },

    #' @description Standard deviation of each intervention between rounds
    sigma = function() {
      if (is.null(private$.sigma)) {
        rho <- private$rho()
        private$.sigma <- sqrt(rho / (1 - rho))
        names(private$.sigma) <- private$interventions
      }
      private$.sigma
    },

    #' @description multivariate norm draws for these parameters
    mvnorm = function() {
      if (is.null(private$.mvnorm)) {
        rho <- private$rho()
        sigma <- self$sigma()
        V <- outer(sigma, sigma) * rho
        diag(V) <- sigma ^ 2
        private$.mvnorm <- MASS::mvrnorm(
          private$population,
          rep(0, length(private$interventions)),
          V
        )
        dimnames(private$.mvnorm)[[2]] <- private$interventions
      }
      private$.mvnorm
    }
  )
)

get_correlation_parameters <- function(parameters) {
  CorrelationParameters$new(parameters)
}

#' @title Sample a population to intervene in given the correlation parameters
#' @param target a vector of individual indices to sample from
#' @param intervention name of the intervention
#' @param p the probability of being selected
#' @param c_param correlation parameters
#' @importFrom stats qnorm
sample_intervention <- function(target, intervention, p, c_param) {
  sigma_squared <- c_param$sigma()[[intervention]]^2
  sd <- sqrt(1 + sigma_squared)
  u0 <- -qnorm(p, 0) * sd
  z <- rnorm(length(target))
  u0 + c_param$mvnorm()[target, intervention] + z < 0
}
