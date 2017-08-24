#' Prepare RNA-Seq RMarkdown Template
#'
#' @rdname prepareRNASeqTemplate
#' @name prepareRNASeqTemplate
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
    prepareTemplate(dir = file.path(url, "downloads"))
})
