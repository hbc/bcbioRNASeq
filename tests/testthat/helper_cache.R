invisible(lapply(
    X = c(
        "bcb_invalid.rda"
    ),
    FUN = function(file, url) {
        if (!file.exists(file)) {
            utils::download.file(
                url = paste(url, file, sep = "/"),
                destfile = file
            )
        }
    },
    url = bcbioRNASeqCacheURL
))
