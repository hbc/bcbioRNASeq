cacheURL <- "http://bcbiornaseq.seq.cloud"
files <- c(
    "sample_metadata.csv"
)
mapply(
    FUN = function(cacheURL, file, envir) {
        # Download file to testthat directory
        if (!file.exists(file)) {
            download.file(
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
