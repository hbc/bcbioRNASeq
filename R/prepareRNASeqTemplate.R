# FIXME This fails during `devtools::run_examples()` because of package path.
# Error in system.file("rmarkdown/shared", package = package, mustWork = TRUE) :
#   no file found



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
