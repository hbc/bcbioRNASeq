#' Prepare RNA-Seq R Markdown Template
#'
#' @name prepareRNASeqTemplate
#' @author Michael Steinbaugh
#' @inherit basejump::prepareTemplate
#' @export
#'
#' @examples
#' x <- prepareRNASeqTemplate()
#' x
#'
#' ## Clean up.
#' unlink(c(
#'     "_footer.Rmd",
#'     "_header.Rmd",
#'     "_output.yaml",
#'     "_setup.R",
#'     "bibliography.bib"
#' ))
prepareRNASeqTemplate <- function(overwrite = FALSE) {
    prepareTemplate(package = "bcbioRNASeq", overwrite = overwrite)
}
