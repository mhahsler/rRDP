.onLoad <- function(libname, pkgname) {
    .jpackage(pkgname, morePaths = system.file("java/rdp_classifier.jar", package = "rRDP"), lib.loc = libname)
}
