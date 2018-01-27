#' Prepare RNA-Seq R Markdown Template
#'
#' @rdname prepareRNASeqTemplate
#' @name prepareRNASeqTemplate
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return No value.
#'
#' @examples
#' \dontrun{
#' prepareRNASeqTemplate()
#' }
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
