#' @title Preset parameters for the An. gambiae s.s vector
#' @details Default parameters:
#' species: "gamb"
#' blood_meal_rates: 0.3333333
#' Q0: 0.92
#' endophily: 0.813
#' rn: 0.56
#' rnm: 0.24
#' dn0: 0.533
#' phi_bednets: 0.89
#' rs: 0.2
#' phi_indoors: 0.97
#' mum: 0.132"
#' @export
gamb_params <- list(
  species = 'gamb',
  blood_meal_rates = 1/3,
  Q0 = .92,
  endophily = .813,
  rn = .56,
  rnm = .24,
  dn0 = .533,
  phi_bednets = .89,
  rs = .2,
  phi_indoors = .97,
  mum = .132
)

#' @title Preset parameters for the An. arabiensis vector
#' @details Default parameters:
#' species: "arab"
#' blood_meal_rates: 0.3333333
#' Q0: 0.71
#' endophily: 0.422
#' rn: 0.46
#' rnm: 0.1
#' dn0: 0.533
#' phi_bednets: 0.9
#' rs: 0.2
#' phi_indoors: 0.96
#' mum: 0.132
#' @export
arab_params <- list(
  species = 'arab',
  blood_meal_rates = 1/3,
  Q0 = .71,
  endophily = .422,
  rn = .46,
  rnm = .1,
  dn0 = .533,
  phi_bednets = .9,
  rs = .2,
  phi_indoors = .96,
  mum = .132
)

#' @title Preset parameters for the An. funestus vector
#' @details Default parameters:
#' species: "fun"
#' blood_meal_rates: 0.3333333
#' Q0: 0.94
#' endophily: 0.813
#' rn: 0.56
#' rnm: 0.1
#' dn0: 0.533
#' phi_bednets: 0.9
#' rs: 0.2
#' phi_indoors: 0.98
#' mum: 0.112
#' @export
fun_params <- list(
  species = 'fun',
  blood_meal_rates = 1/3,
  Q0 = .94,
  endophily = .813,
  rn = .56,
  rnm = .1,
  dn0 = .533,
  phi_bednets = .9,
  rs = .2,
  phi_indoors = .98,
  mum = .112
)

#' @title Parameterise the mosquito species to use in the model
#'
#' @param parameters the model parameters
#' @param species a list of species presets
#' @param proportions a vector of proportions for each species
#' @export
set_species <- function(parameters, species, proportions) {
  if (length(species) != length(proportions)) {
    stop('You must give proportions for each species')
  }
  if (!approx_sum(proportions, 1)) {
    stop('Proportions do not sum to 1')
  }
  keys <- c(
    'species',
    'blood_meal_rates',
    'Q0',
    'endophily',
    'rn',
    'rnm',
    'dn0',
    'phi_bednets',
    'rs',
    'phi_indoors',
    'mum'
  )
  for (key in keys) {
    parameters[[key]] <- rep(NA, length(species))
  }
  for (v in seq_along(species)) {
    for (key in keys) {
      parameters[[key]][[v]] <- species[[v]][[key]]
    }
  }
  parameters$species_proportions <- proportions
  parameters
}
