#' bcbioRNADataSet Object Class (Legacy)
#'
#' This class has been superceded by [bcbioRNASeq] and will be formally
#' deprecated in a future release.
#'
#' @author Lorena Pantano, Michael Steinbaugh
#' @keywords internal
#'
#' @slot callers [SimpleList] containing additional bcbio run data with
#'   dimensions that don't match the count matrix.
#'
#' @export
bcbioRNADataSet <- setClass(
    "bcbioRNADataSet",
    contains = "SummarizedExperiment",
    slots = c(callers = "SimpleList")
)
setValidity("bcbioRNADataSet", function(object) TRUE)
