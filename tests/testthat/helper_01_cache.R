cacheURL <- "http://bcbiornaseq.seq.cloud"
files <- c(
    "bcb_invalid.rda"
)
mapply(
    FUN = function(cacheURL, file, envir) {
        if (!file.exists(file)) {
            utils::download.file(
                url = paste(cacheURL, file, sep = "/"),
                destfile = file)
        }
        # Load R Data file
        if (grepl("\\.rda$", file)) {
            message(paste("Loading", file))
            load(file, envir = envir)
        }
    },
    file = files,
    MoreArgs = list(cacheURL = cacheURL, envir = environment())
)
