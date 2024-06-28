## code to prepare parasite parameters
parasite_csv <- read.csv("data-raw/parasite_parameters.csv")
parasite_parameters <- lapply(unique(parasite_csv$parasite), function(x){
  parameters <- as.list(parasite_csv[parasite_csv$parasite==x,c("default")])
  names(parameters) <- parasite_csv[parasite_csv$parasite==x,c("parameter")]
  return(parameters)
})
names(parasite_parameters) <- unique(parasite_csv$parasite)
usethis::use_data(parasite_parameters, overwrite = TRUE)