if (!isTRUE(hasInternet())) {
    warning("No Internet connection detected.")
    return()
}
dir.create("cache", showWarnings = FALSE)
files <- "bcb_invalid.rda"
mapply(
    FUN = function(remoteDir, file, envir) {
        destfile <- file.path("cache", file)
        if (!file.exists(destfile)) {
            utils::download.file(
                url = paste(remoteDir, file, sep = "/"),
                destfile = destfile
            )
        }
    },
    file = files,
    MoreArgs = list(
        remoteDir = bcbioRNASeqTestsURL,
        envir = environment()
    )
)
