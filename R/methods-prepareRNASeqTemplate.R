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
#'     "bibliography.bib",
#'     "setup.R"
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
            sourceDir = system.file("rmarkdown/shared",
                                    package = "bcbioRNASeq"))
    })
