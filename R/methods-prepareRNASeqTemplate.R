#' Prepare RNA-Seq R Markdown Template
#'
#' @rdname prepareRNASeqTemplate
#' @name prepareRNASeqTemplate
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return No value.
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
NULL


# Methods ======================================================================
#' @rdname prepareRNASeqTemplate
#' @importFrom bcbioBase prepareTemplate
#' @export
setMethod(
    "prepareRNASeqTemplate",
    signature("missing"),
    function(object) {
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
                package = "bcbioRNASeq")
        )
    })
