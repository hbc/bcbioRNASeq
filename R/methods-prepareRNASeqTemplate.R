#' Prepare RNA-Seq RMarkdown Template
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



# Methods ====
#' @rdname prepareRNASeqTemplate
#' @export
setMethod("prepareRNASeqTemplate", "missing", function(object) {
    prepareTemplate(
        sourceDir = system.file("rmarkdown/shared", package = "bcbioRNASeq"))
})
