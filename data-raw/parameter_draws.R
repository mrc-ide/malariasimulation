## code to prepare parameter draws

pd <- read.csv("data-raw/parameter_draws.csv")

# Split individual draws into a list of named parameter lists
draws <- 1:max(pd$draw)
parameter_draws <- list()
for(draw in draws){
  pars <- as.list(pd[pd$draw == draw, "value"])
  names(pars) <- pd[pd$draw == draw, "parameter"]
  parameter_draws[[draw]] <- pars
}
names(parameter_draws) <- draws

usethis::use_data(parameter_draws, overwrite = TRUE)
