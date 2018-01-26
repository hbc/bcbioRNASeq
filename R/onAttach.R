.onAttach <- function(libname, pkgname) {
    packages <- c(
        "bcbioBase",
        "SummarizedExperiment",
        "viridis",
        "DESeq2",
        "DEGreport"
    )
    lapply(packages, function(package) {
        if (!package %in% (.packages())) {
            attachNamespace(package)
        }
    })
    invisible()
}
