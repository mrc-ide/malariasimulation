#' @title Parameterise a forced (time-varying) EIR
#'
#' @description When a forced EIR is set, the mosquito component of the model is
#' not simulated. Instead, the number of infectious bites delivered to humans at
#' each timestep is driven directly from the user-supplied EIR. The overdispersed
#' distribution of bites among individuals (heterogeneous biting by age and
#' individual biting rate) is preserved.
#'
#' The EIR is specified as a step function: each value applies from its
#' corresponding timestep until the next change point (and is 0 before the first
#' change point, so set the first timestep to 1 to force the EIR from the start).
#' EIR is given in units of infectious bites per adult person per year - the same
#' per-adult convention as the default (`EIR_population_input = "adult"`)
#' `init_EIR` argument of \code{\link{set_equilibrium}}. Forcing `EIR = X`
#' therefore reproduces a `set_equilibrium(init_EIR = X)` run. The total EIR is
#' split across mosquito species according to `species_proportions` (i.e. by
#' biting share, not by any per-species difference in transmission efficiency).
#'
#' Because the mosquito population is not simulated, the prescribed EIR is applied
#' directly and does not respond to processes that would ordinarily feed back
#' through the vector population, so \code{set_forced_eir} emits a warning. In
#' particular, when combining a forced EIR with interventions:
#' \itemize{
#'   \item Interventions that act \strong{directly on humans} (drugs, MDA/SMC/PMC,
#'   vaccines) work as usual.
#'   \item \strong{Personal-protection} effects of vector control
#'   (\code{\link{set_bednets}}, \code{\link{set_spraying}}) still
#'   \emph{redistribute} the fixed number of infectious bites away from protected
#'   individuals, but do not change the total EIR. If your prescribed EIR was
#'   estimated from a setting that already had these interventions, applying them
#'   again would double-count their effect - supply the pre-intervention EIR.
#'   \item Interventions that act \strong{only through the mosquito population}
#'   (\code{\link{set_carrying_capacity}}, larval source management) have
#'   \strong{no effect} under a forced EIR.
#' }
#' There is also no entomological lag: the forced EIR is the infectious-bite rate
#' experienced at timestep `t`, so responses to step changes appear about
#' `parameters$de` days earlier than in a mosquito-simulated run. Seasonality is
#' likewise bypassed and must be supplied directly as a fine-grained EIR series.
#' Human-side dynamics (immunity, superinfection, and P. vivax relapse from the
#' hypnozoite reservoir) continue to operate normally.
#'
#' \code{\link{set_equilibrium}} should still be used to initialise the human
#' population. A forced EIR is not supported in metapopulation simulations, and
#' \code{force_EIR} must be kept consistent if a run is resumed from a saved state.
#'
#' @param parameters the model parameters
#' @param timesteps vector of timesteps at which the EIR changes
#' @param eir vector of EIR values (infectious bites per adult person per year),
#' the same length as `timesteps`
#'
#' @seealso \code{\link{set_equilibrium}}
#' @export
set_forced_eir <- function(parameters, timesteps, eir) {
  stopifnot(length(eir) == length(timesteps))
  stopifnot(all(eir >= 0))
  stopifnot(min(timesteps) > 0)
  stopifnot(!is.unsorted(timesteps))

  warning(
    paste(
      "Forcing the EIR bypasses the simulated mosquito population: the",
      "prescribed EIR is applied directly and will not respond to interventions",
      "or other processes that would ordinarily feed back through the vector",
      "population. Ensure the forced EIR already reflects the transmission you",
      "intend to model. See ?set_forced_eir for how interventions interact with",
      "a forced EIR."
    ),
    call. = FALSE
  )

  if (min(timesteps) > 1) {
    warning(
      paste0(
        "The first forced EIR change point is at timestep ", min(timesteps),
        "; the EIR is 0 before then. Set the first timestep to 1 to force the ",
        "EIR from the start of the simulation."
      ),
      call. = FALSE
    )
  }

  parameters$force_EIR <- TRUE
  parameters$force_EIR_timesteps <- timesteps
  parameters$force_EIR_values <- eir
  parameters
}

#' @title Look up the forced EIR at a given timestep
#'
#' @description Returns the total forced EIR (infectious bites per person per
#' year) applicable at `timestep`, as a piecewise-constant step function of the
#' timesteps and values set by \code{\link{set_forced_eir}}. Returns 0 before the
#' first change point.
#'
#' @param parameters the model parameters
#' @param timestep the current timestep
#' @noRd
get_forced_eir <- function(parameters, timestep) {
  previous <- which(parameters$force_EIR_timesteps <= timestep)
  if (length(previous) == 0) {
    return(0)
  }
  parameters$force_EIR_values[[max(previous)]]
}

#' @title Convert a forced EIR into the native per-timestep EIR driver
#'
#' @description Converts the user-supplied EIR (infectious bites per adult person
#' per year) for a species into the native per-timestep EIR value used to drive
#' biting in `simulate_bites`.
#'
#' The input EIR follows the same per-adult convention as the `init_EIR` argument
#' of \code{\link{set_equilibrium}} (its default `EIR_population_input = "adult"`).
#' We first convert it to a whole-population-average per-person rate by dividing by
#' the same adult scaling factor `set_equilibrium` uses
#' (\code{calculate_population_to_adult_EIR_scalar}), then scale to a
#' population-total daily driver. The conversion is chosen so that the downstream
#' `expected_bites <- species_eir * mean(psi)` calculation yields the intended
#' population-total number of infectious bites per timestep, and so that forcing
#' `EIR = X` reproduces a `set_equilibrium(init_EIR = X)` run.
#'
#' @param parameters the model parameters
#' @param species the species index
#' @param psi the age-based relative biting rate for each human
#' @param timestep the current timestep
#' @noRd
forced_species_eir <- function(parameters, species, psi, timestep) {
  adult_eir_pppy <- get_forced_eir(parameters, timestep) *
    parameters$species_proportions[[species]]
  # Convert per-adult EIR to a whole-population-average per-person rate using the
  # same scaling set_equilibrium applies (adult = population * omega_age).
  omega_age <- calculate_population_to_adult_EIR_scalar(parameters, adult_eir_pppy)
  population_eir_pppy <- adult_eir_pppy / omega_age
  population_eir_pppy / 365 * parameters$human_population / mean(psi)
}
