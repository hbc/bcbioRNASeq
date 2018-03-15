#' Prepare RNA-Seq R Markdown Template
#'
#' @name prepareRNASeqTemplate
#' @family R Markdown Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase prepareTemplate
#'
#' @inheritParams general
#'
#' @return No value.
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
prepareRNASeqTemplate <- function() {
    prepareTemplate(
        c(
            "_footer.Rmd",
            "_header.Rmd",
            "_output.yaml",
            "_setup.R",
            "bibliography.bib"
        ),
        sourceDir = system.file(
            "rmarkdown/shared",
            package = "bcbioRNASeq"
        )
    )
}
