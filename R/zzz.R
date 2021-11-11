
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(c(strwrap(
      paste("NOTE: tergmLite functionality has largely been wrapped into EpiModel
            with EpiModel >= v2.2.0. Some updates may be needed for EpiModel extension
            models or packages that modify the network resimulation modules in EpiModel.
            File a help issue at https://github.com/EpiModel/EpiModel for assistance.",
            sep = "")), ""), collapse = "\n"))
}
