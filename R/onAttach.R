.onAttach <- function(libname, pkgname) {
    imports <- c(
        "bcbioBase",
        "SummarizedExperiment",
        "viridis",
        "DESeq2",
        "DEGreport"
    )
    invisible(lapply(
        X = imports,
        FUN = require,
        character.only = TRUE
    ))
}
