#' Download Dependency File
#'
#' If the required dependency file isn't present, download latest version from
#' the [HBC website](http://bioinformatics.sph.harvard.edu/).
#'
#' @author Michael Steinbaugh
#'
#' @param file *Optional*. File name. If `NULL` (default), download the default
#'   dependency files for a new experiment.
#'
#' @export
download <- function(file = NULL) {
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
        dl(c("_footer.Rmd",
             "_header.Rmd",
             "_output.yaml",
             "bcbioRnaseq.bib",
             "load.R",
             "setup.R"))
    }
}
