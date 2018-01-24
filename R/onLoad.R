.onLoad <- function(libname, pkgname) {
    packages <-
        c("bcbioBase",
          "SummarizedExperiment",
          "viridis",
          "DESeq2",
          "DEGreport")
    lapply(seq_along(packages), function(a) {
        if (!packages[[a]] %in% (.packages())) {
            attachNamespace(packages[[a]])
        }
    })
    invisible()
}
