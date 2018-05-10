#' Prepare RNA-Seq R Markdown Template
#'
#' @name prepareRNASeqTemplate
#' @family R Markdown Functions
#' @author Michael Steinbaugh
#'
#' @inherit bcbioBase::prepareTemplate
#'
#' @export
#'
#' @examples
#' prepareRNASeqTemplate()
#'
#' # Clean up
#' unlink(c(
#'     "_footer.Rmd",
#'     "_header.Rmd",
#'     "_output.yaml",
#'     "_setup.R",
#'     "bibliography.bib"
#' ))
prepareRNASeqTemplate <- function(overwrite = FALSE) {
    prepareTemplate(
        file = c(
            "_footer.Rmd",
            "_header.Rmd",
            "_output.yaml",
            "_setup.R",
            "bibliography.bib"
        ),
        sourceDir = system.file(
            "rmarkdown/shared",
            package = "bcbioRNASeq"
        ),
        overwrite = overwrite
    )
}
