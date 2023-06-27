## code to prepare parasite parameters
parasite_parameters <- read.csv("data-raw/parasite_parameters.csv")
usethis::use_data(parasite_parameters, overwrite = TRUE)
