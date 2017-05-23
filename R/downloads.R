#' File downloads
#'
#' If file isn't present, download latest version from the
#' [HBC website](http://bioinformatics.sph.harvard.edu/).
#'
#' File download utility functions for RMarkdown knit reports.
#'
#' @author Michael Steinbaugh
#'
#' @param file File name. If `NULL` (default), download the default dependency
#'   files for a consult.
#'
#' @export
downloads <- function(file = NULL) {
    dl <- function(file) {
        sapply(seq_along(file), function(a) {
            if (!file.exists(file[a])) {
                download.file(
                    file.path("http://bioinformatics.sph.harvard.edu",
                              "bcbioRnaseq", "downloads", file[a]),
                    destfile = file[a])
            }
        }) %>% invisible
    }

    if (!is.null(file)) {
        dl(file)
    } else {
        # HBC project defaults
        dl(c("_output.yaml",
             "bcbioRnaseq.bib",
             "footer.Rmd",
             "header.Rmd",
             "setup.R"))
    }
}
