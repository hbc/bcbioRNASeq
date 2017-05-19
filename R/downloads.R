#' File downloads
#'
#' If file isn't present, download latest version from the
#' [HBC website](http://bioinformatics.sph.harvard.edu/).
#'
#' File download utility functions for RMarkdown knit reports.
#'
#' @rdname downloads
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param file File name.



#' @rdname downloads
hbc_download <- function(file) {
    sapply(seq_along(file), function(a) {
        if (!file.exists(file[a])) {
            download.file(
                file.path("http://bioinformatics.sph.harvard.edu",
                          "bcbioRnaseq", "downloads", file[a]),
                destfile = file[a])
        }
    }) %>% invisible
}



#' @rdname downloads
download_rmarkdown_files <- function() {
    hbc_download(
        c("_output.yaml",
          "bcbioRnaseq.bib",
          "setup.R"))
}
