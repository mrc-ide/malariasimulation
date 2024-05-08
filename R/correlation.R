INTS <- c(
  'pmc',
  'pev',
  'mda',
  'smc',
  'tbv',
  'bednets',
  'spraying'
)

#' Class: Correlation parameters
#'
#' This class implements functionality that allows interventions to be
#' correlated, positively or negatively. By default, interventions are applied
#' independently and an individual's probability of receiving two interventions
#' (either two separate interventions or two rounds of the same one) is the
#' product of the probability of receiving each one.
#'
#' By setting a positive correlation between two interventions, we can make it
#' so that the individuals that receive intervention A are more likely to
#' receive intervention B. Conversely, a negative correlation will make it such
#' that individuals that receive intervention A are less likely to also receive
#' intervention B.
#'
#' Broadly speaking, the implementation works by assigning at startup a weight
#' to each individual and intervention pair, reflecting how likely an individual
#' is to receive that intervention. Those weights are derived stochastically
#' from the configured correlation parameters.
#'
#' For a detailed breakdown of the calculations, see Protocol S2 of
#' Griffin et al. (2010).
CorrelationParameters <- R6::R6Class(
  'CorrelationParameters',
  private = list(
    interventions = list(),
    n_ints = function() length(private$interventions),
    population = 0,
    rho_matrix = NULL,
    rho = function() diag(private$rho_matrix),
    .sigma = NULL,
    .mvnorm = NULL,

    #' Derive the mvnorm from the configured correlations.
    #'
    #' If a \code{restored_mvnorm} is specified, its columns (corresponding to
    #' restored interventions) will be re-used as is. Missing columns (for new
    #' interventions) are derived in accordance with the restored data.
    calculate_mvnorm = function(restored_mvnorm = matrix(ncol=0, nrow=private$population)) {
      sigma <- self$sigma()
      V <- outer(sigma, sigma) * private$rho_matrix
      diag(V) <- sigma ^ 2

      restored_interventions <- match(colnames(restored_mvnorm), private$interventions)
      new_interventions <- setdiff(seq_along(private$interventions), restored_interventions)

      mvnorm <- matrix(
        nrow = private$population,
        ncol = length(private$interventions),
        dimnames = list(NULL, private$interventions)
      )
      mvnorm[,restored_interventions] <- restored_mvnorm
      if (length(new_interventions) > 0) {
        mvnorm[,new_interventions] <- rcondmvnorm(
          private$population,
          mean = rep(0, length(private$interventions)),
          sigma = V,
          given = restored_mvnorm,
          dependent.ind = new_interventions,
          given.ind = restored_interventions
        )
      }

      mvnorm
    }
  ),
  public = list(

    #' @description initialise correlation parameters
    #' @param population popularion size
    #' @param interventions character vector with the name of enabled interventions
    initialize = function(population, interventions) {
      private$population <- population
      private$interventions <- interventions

      # Initialise a rho matrix for our interventions
      n_ints <- private$n_ints()
      private$rho_matrix <- matrix(
        0,
        nrow = n_ints,
        ncol = n_ints,
        dimnames = list(private$interventions, private$interventions)
      )
    },

    #' @description Add rho between rounds
    #' @param int string representing the intervention to update
    #' @param rho value between 0 and 1 representing the correlation between rounds of
    #' the intervention
    inter_round_rho = function(int, rho) {
      stopifnot(is.null(private$.sigma) && is.null(private$.mvnorm))

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
    },

    #' @description Add rho between interventions
    #' @param int_1 string representing the first intervention
    #' @param int_2 string representing the second intervention (intechangable
    #' with int_1)
    #' @param rho value between -1 and 1 representing the correlation between rounds of
    #' the intervention
    inter_intervention_rho = function(int_1, int_2, rho) {
      stopifnot(is.null(private$.sigma) && is.null(private$.mvnorm))

      if (!(int_1 %in% private$interventions)) {
        stop(paste0('invalid intervention name: ', int_1))
      }
      if (!(int_2 %in% private$interventions)) {
        stop(paste0('invalid intervention name: ', int_2))
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
        private$.mvnorm <- private$calculate_mvnorm()
      }
      private$.mvnorm
    },

    #' @description Save the correlation state.
    save_state = function() {
      # mvnorm is sampled at random lazily on its first use. We need to save it
      # in order to restore the same value when resuming the simulation,
      # otherwise we would be drawing a new, probably different, value.
      # The rest of the object is derived deterministically from the parameters
      # and does not need saving.
      list(mvnorm=self$mvnorm())
    },

    #' @description Restore the correlation state.
    #'
    #' Only the randomly drawn weights are restored. The object needs to be
    #' initialized with the same rhos.
    #'
    #' @param timestep the timestep at which simulation is resumed. This
    #' parameter's value is ignored, it only exists to conform to a uniform
    #' interface.
    #' @param state a previously saved correlation state, as returned by the
    #' save_state method.
    restore_state = function(timestep, state) {
      stopifnot(is.null(private$.sigma) && is.null(private$.mvnorm))
      private$.mvnorm <- private$calculate_mvnorm(state$mvnorm)
    }
  )
)

#' @title Get default correlation parameters
#' @description returns a `CorrelationParameters` object for you edit. By
#' default, all correlations are set to 0
#'
#' @param parameters model parameters
#' @export
#' @examples
#' 
#' # get the default model parameters
#' parameters <- get_parameters()
#' 
#' # Set some vaccination strategy
#' parameters <- set_mass_pev(
#'   parameters,
#'   profile = rtss_profile,
#'   timesteps = 100,
#'   coverages = .9,
#'   min_wait = 0,
#'   min_ages = 100,
#'   max_ages = 1000,
#'   booster_spacing = numeric(0),
#'   booster_coverage = numeric(0),
#'   booster_profile = NULL
#' )
#' 
#' # Set some smc strategy
#' parameters <- set_drugs(parameters, list(SP_AQ_params))
#' parameters <- set_smc(
#'   parameters,
#'   drug = 1,
#'   timesteps = 100,
#'   coverages = .9,
#'   min_age = 100,
#'   max_age = 1000
#' )
#' 
#' # Correlate the vaccination and smc targets
#' correlations <- get_correlation_parameters(parameters)
#' correlations$inter_intervention_rho('pev', 'smc', 1)
#' 
#' # Correlate the rounds of smc
#' correlations$inter_round_rho('smc', 1)
#' 
#' # You can now pass the correlation parameters to the run_simulation function
get_correlation_parameters <- function(parameters) {
  # Find a list of enabled interventions
  enabled <- vlapply(INTS, function(name) parameters[[name]])

  CorrelationParameters$new(parameters$human_population, INTS[enabled])
}

#' @title Sample a population to intervene in given the correlation parameters
#' @param target a vector of individual indices to sample from
#' @param intervention name of the intervention
#' @param p the probability of being selected
#' @param correlations correlation parameters
#' @importFrom stats qnorm
#' @noRd
sample_intervention <- function(target, intervention, p, correlations) {
  sigma_squared <- correlations$sigma()[[intervention]]^2
  sd <- sqrt(1 + sigma_squared)
  u0 <- -qnorm(p, 0) * sd
  z <- rnorm(length(target))
  u0 + correlations$mvnorm()[target, intervention] + z < 0
}

#' Simulate from a conditional multivariate normal distribution.
#'
#' Given a multidimensional variable Z which follows a multivariate normal
#' distribution, this function allows one to draw samples for a subset of Z,
#' while putting conditions on the values of the rest of Z.
#'
#' This effectively allows one to grow a MVN distributed matrix (with columns as
#' the dimensions and a row per sampled vector), adding new dimensions after the
#' fact. The existing columns are used as the condition set on the distribution,
#' and the values returned by this function are used as the new dimensions.
#'
#' The maths behind the implementation are described in various online sources:
#' - https://statproofbook.github.io/P/mvn-cond.html
#' - https://www.stats.ox.ac.uk/~doucet/doucet_simulationconditionalgaussian.pdf
#' - https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
#'
#' @param n the number of samples to simulate
#' @param mean the mean vector of the distribution, including both given and
#' dependent variables
#' @param sigma the variance-covariance matrix of the distribution, including
#' both given and dependent variables
#' @param given a matrix of given values used as conditions when simulating the
#' distribution. The matrix should have \code{n} rows, each one specifying a
#' different set of values for the given variables.
#' @param dependent.ind the indices within \code{mean} and \code{sigma} of the
#' variables to simulate.
#' @param given.ind the indices within \code{mean} and \code{sigma} of the
#' variables for which conditions are given. The length of this vector must be
#' equal to the number of columns of the \code{given} matrix. If empty or NULL,
#' this function is equivalent to simulating from an unconditional multivariate
#' normal distribution.
#' @return a matrix with \code{n} rows and \code{length(dependent.ind)} columns,
#' containing the simulated value.
#' @importFrom MASS mvrnorm
#' @noRd
rcondmvnorm <- function(n, mean, sigma, given, dependent.ind, given.ind) {
  stopifnot(length(mean) == nrow(sigma))
  stopifnot(length(mean) == ncol(sigma))
  stopifnot(nrow(given) == n)
  stopifnot(ncol(given) == length(given.ind))

  sigma11 <- sigma[dependent.ind, dependent.ind, drop=FALSE]
  sigma12 <- sigma[dependent.ind, given.ind, drop=FALSE]
  sigma21 <- sigma[given.ind, dependent.ind, drop=FALSE]
  sigma22 <- sigma[given.ind, given.ind, drop=FALSE]

  if (all(sigma22 == 0)) {
    # This covers two cases: there were no given variables and therefore their
    # variance-covariance matrix is empty, or there were given variables but
    # they had a variance of zero. The general formula can't support the latter
    # case since it tries to invert the matrix, but we can safely ignore the
    # values since they are all equal to their mean and don't influence the
    # dependent variables.
    #
    # In both cases we revert to a standard MVN with no condition.
    mvrnorm(n, mean[dependent.ind], sigma11)
  } else {
    # Available implementations of the conditional multivariate normal assume
    # every sample is drawn using the same condition on the given variables.
    # This is not true in our usecase, where every individual has already had an
    # independent vector of values drawn for the given variable. We are
    # effectively drawing from as many different distributions as there are
    # individuals. Thankfully the same conditional covariance matrix can be
    # used for all the distributions, only the mean vector needs to be
    # different. We draw the underlying samples from the MVN at mean 0, and
    # offset that later on a per-individual basis.
    #
    # To work over all the vectors directly they need to be as columns, which
    # is why we start by transposing `given`. R will recycle the `m` matrix and
    # `mean` vectors across all the columns. The last step is to transpose the
    # result back into the expected configuration.

    m <- sigma12 %*% solve(sigma22)
    residual <- t(given) - mean[given.ind]
    cond_mu <- t(m %*% residual + mean[dependent.ind])
    cond_sigma <- sigma11 - m %*% sigma21

    samples <- mvrnorm(n, rep(0, length(dependent.ind)), cond_sigma)
    samples + cond_mu
  }
}

used_intervention <- function(variable, timestep, window) {
  variable$get_index_of(set=-1)$not()$and(
    variable$get_index_of(a=timestep - window, b=timestep)
  )
}

create_combined_intervention_rendering_process <- function(
  int_1,
  variable_1,
  int_2,
  variable_2,
  window,
  renderer
) {
  name <- paste0('n_combined_', int_1, '_', int_2)
  renderer$set_default(name, 0)
  function (timestep) {
    n <- used_intervention(variable_1, timestep, window)$and(
      used_intervention(variable_2, timestep, window)
    )$size()
    renderer$render(
      name,
      n,
      timestep
    )
  }
}
