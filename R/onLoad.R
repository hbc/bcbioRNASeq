.onLoad <- function(libname, pkgname) {
    if (!"annotables" %in% (.packages())) {
        attachNamespace("annotables")
    }
    invisible()
}
