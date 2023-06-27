
# read.tcsv <- function(file, header=TRUE, sep=",", ...) {
#   
#   n = max(count.fields(file, sep=sep), na.rm=TRUE)
#   x = readLines(file)
#   
#   .splitvar = function(x, sep, n) {
#     var = unlist(strsplit(x, split=sep))
#     length(var) = n
#     return(var)
#   }
#   
#   x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
#   x = apply(x, 1, paste, collapse=sep) 
#   out = read.csv(text=x, sep=sep, header=header, ...)
#   return(out)
#   
# }
# 
# 
# parasite_parameters <- read.tcsv("data-raw/parasite_parameters_long_form.csv", nrows = 3)
# usethis::use_data(parasite_parameters, overwrite = TRUE)

parasite_parameters <- read.csv("data-raw/parasite_parameters_long_form_references.csv")
usethis::use_data(parasite_parameters, overwrite = TRUE)

## code to prepare parasite parameters
# parasite_parameters <- read.csv("data-raw/parasite_parameters.csv", nrows = 3)
# usethis::use_data(parasite_parameters, overwrite = TRUE)

