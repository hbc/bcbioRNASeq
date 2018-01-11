.onLoad <- function(libname, pkgname) {
    packages <-
        c("SummarizedExperiment",
          "bcbioBase",
          "DESeq2",
          "DEGreport")
    lapply(seq_along(packages), function(a) {
        if (!packages[[a]] %in% (.packages())) {
            attachNamespace(packages[[a]])
        }
    })
    invisible()
}
