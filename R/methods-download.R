#' Download Dependency File
#'
#' If the required dependency file isn't present, download latest version from
#' the [HBC website](http://bioinformatics.sph.harvard.edu/).
#'
#' File download utility for RMarkdown knit reports.
#'
#' @rdname download
#' @name download
#'
#' @param object *Optional*. File name. If `NULL` (default), download the
#'   default dependency files for a new experiment.
#'
#' @examples
#' download("setup.R")
NULL



# Constructors ====
.download <- function(file) {
    sapply(seq_along(file), function(a) {
        if (!file.exists(file[[a]])) {
            download.file(
                file.path(packageURL, "downloads", file[[a]]),
                destfile = file[[a]])
        }
    }) %>% invisible
}



# Methods ====
#' @rdname download
#' @export
setMethod("download", "missing", function() {
    .download(
        c("_output.yaml",
          "_footer.Rmd",
          "_header.Rmd",
          "bcbioRnaseq.bib",
          "setup.R"))
})



#' @rdname download
#' @export
setMethod("download", "character", function(object) {
    .download(object)
})
