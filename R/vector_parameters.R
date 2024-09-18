#' @title Preset parameters for the An. gambiae s.s vector
#' @details Default parameters:
#' species: "gamb"
#' blood_meal_rates: 0.3333333
#' foraging_time: .69
#' Q0: 0.92
#' phi_bednets: 0.85
#' phi_indoors: 0.90
#' mum: 0.132
#'
#' parameters from:
#' https://www.pnas.org/content/pnas/early/2019/07/02/1820646116.full.pdf
#' @export
gamb_params <- list(
  species = 'gamb',
  blood_meal_rates = 1/3,
  foraging_time = .69,
  Q0 = .92,
  phi_bednets = .85,
  phi_indoors = .90,
  mum = .132
)

#' @title Preset parameters for the An. arabiensis vector
#' @details Default parameters:
#' species: "arab"
#' blood_meal_rates: 0.3333333
#' foraging_time: .69
#' Q0: 0.71
#' phi_bednets: 0.8
#' phi_indoors: 0.86
#' mum: 0.132
#'
#' parameters from:
#' https://www.pnas.org/content/pnas/early/2019/07/02/1820646116.full.pdf
#' @export
arab_params <- list(
  species = 'arab',
  blood_meal_rates = 1/3,
  foraging_time = .69,
  Q0 = .71,
  phi_bednets = .8,
  phi_indoors = .86,
  mum = .132
)

#' @title Preset parameters for the An. funestus vector
#' @details Default parameters:
#' species: "fun"
#' blood_meal_rates: 0.3333333
#' foraging_time: .69
#' Q0: 0.94
#' phi_bednets: 0.78
#' phi_indoors: 0.87
#' mum: 0.112
#'
#' parameters from:
#' https://www.pnas.org/content/pnas/early/2019/07/02/1820646116.full.pdf
#' @export
fun_params <- list(
  species = 'fun',
  blood_meal_rates = 1/3,
  foraging_time = .69,
  Q0 = .94,
  phi_bednets = .78,
  phi_indoors = .87,
  mum = .112
)

#' @title Preset parameters for the An. stephensi vector
#' @details Default parameters:
#' species: "steph"
#' blood_meal_rates: 0.3333333
#' foraging_time: .69
#' Q0: 0.21
#' phi_bednets: 0.57
#' phi_indoors: 0.37
#' mum: 0.112
#'
#' parameters reference:
#' https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-022-02324-1
#' Values for phi are from:
#' https://github.com/arranhamlet/stephensi_ETH_publication/blob/297352e244f8ed658e8bc3f32be42f011269c7f0/R/functions/cluster_hypercube_sampling_djibouti.R#L14-L21
#' values for Q0: are the average from:
#' ttps://github.com/cwhittaker1000/stephenseasonality/blob/main/data/bionomic_species_all_LHC_100.csv
#' @export
steph_params <- list(
  species = 'steph',
  blood_meal_rates = 1/3,
  foraging_time = .69,
  Q0 = 0.21,
  phi_bednets = 0.52186,
  phi_indoors = 0.4776,
  mum = .112
)

#' @title Preset parameters for the An. koliensis vector
#' @details Default parameters:
#' species: "kol"
#' blood_meal_rates: 0.3333333
#' foraging_time: .68
#' Q0: 0.5
#' phi_bednets: 0.9
#' phi_indoors: 0.5
#' mum: 0.1666667
#'
#' parameters reference:
#' https://github.com/MWhite-InstitutPasteur/Pvivax_IBM/blob/master/koliensis_parameters.txt
#' Value for blood meal rate is assumed in absence of  phi are from:
#' @export
kol_params <- list(
  species = 'kol',
  blood_meal_rates = 1/3,
  foraging_time = 0.68,
  Q0 = 0.5,
  phi_bednets = 0.9,
  phi_indoors = 0.5,
  mum = 0.1666667
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
    'foraging_time',
    'Q0',
    'phi_bednets',
    'phi_indoors',
    'mum'
  )
  for (key in keys) {
    parameters[[key]] <- rep(NA, length(species))
  }
  for (v in seq_along(species)) {
    if (species[[v]]$foraging_time > 1 / species[[v]]$blood_meal_rates) {
      stop(
        sprintf(
          "blood meal time (%f) must be >= foraging time (%f)",
          1 / species[[v]]$blood_meal_rates,
          species[[v]]$foraging_time
        )
      )
    }
    for (key in keys) {
      parameters[[key]][[v]] <- species[[v]][[key]]
    }
  }
  parameters$species_proportions <- proportions
  parameters
}
