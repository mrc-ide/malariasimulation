.onLoad <- function(libname, pkgname){
  headers_version <- individual_headers_version()
  package_version <- as.character(packageVersion("individual"))
  if (!identical(package_version, headers_version)) {
    warning(sprintf(paste(
      "malariasimulation was compiled against version '%s' of the individual,",
      "but version '%s' is currently loaded. This is likely to cause,",
      "malariasimulation to misbehave. You should re-install",
      "malariasimulation. If using devtools, run `devtools::clean_dll` and",
      "then load malariasimulation again."
    ), headers_version, package_version))
  }
}
