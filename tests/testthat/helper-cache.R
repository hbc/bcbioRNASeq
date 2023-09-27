if (!isTRUE(goalie::hasInternet())) {
    warning("No Internet connection detected.")
    return(invisible(NULL))
}
dir.create("cache", showWarnings = FALSE)
files <- c(
    "bcb_fast.rds",
    "fastrnaseq.tar.gz"
)
Map(
    f = function(remoteDir, file, envir) {
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
        "remoteDir" = bcbioRnaSeqTestsUrl,
        "envir" = environment()
    )
)
rm(files)
if (!isTRUE(dir.exists(file.path("cache", "fastrnaseq")))) {
    utils::untar(
        tarfile = file.path("cache", "fastrnaseq.tar.gz"),
        exdir = "cache"
    )
}
stopifnot(dir.exists(file.path("cache", "fastrnaseq")))
