cacheURL <- "http://bcbiornaseq.seq.cloud"
files <- "bcb_invalid.rda"
mapply(
    FUN = function(cacheURL, file, envir) {
        if (!file.exists(file)) {
            utils::download.file(
                url = paste(cacheURL, file, sep = "/"),
                destfile = file)
        }
    },
    file = files,
    MoreArgs = list(cacheURL = cacheURL, envir = environment())
)
