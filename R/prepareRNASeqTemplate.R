#' Prepare RNA-seq R Markdown template
#'
#' @name prepareRNASeqTemplate
#' @author Michael Steinbaugh
#' @inherit basejump::prepareTemplate
#' @export
#'
#' @examples
#' ## x <- prepareRNASeqTemplate()
prepareRNASeqTemplate <- function(overwrite = FALSE) {
    prepareTemplate(package = "bcbioRNASeq", overwrite = overwrite)
}
