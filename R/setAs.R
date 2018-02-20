# TODO Rename to `updateObject()`

#' Coerce Object
#'
#' @rdname coerce
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @param from Class for which the coerce method will perform coercion.
#'
#' @seealso `help(topic = "coerce", package = "methods")`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' se <- as(bcb, "SummarizedExperiment")
#' print(se)
NULL



# Constructors =================================================================
.coerceToSummarizedExperiment <- function(from) {
    to <- new("SummarizedExperiment")
    slot(to, "colData") <- slot(from, "colData")
    slot(to, "assays") <- slot(from, "assays")
    slot(to, "NAMES") <- slot(from, "NAMES")
    slot(to, "elementMetadata") <- slot(from, "elementMetadata")
    slot(to, "metadata") <- slot(from, "metadata")
    validObject(to)
    to
}



# Methods ======================================================================
#' @rdname coerce
#' @name coerce-bcbioRNASeq-SummarizedExperiment
#' @section bcbioRNASeq to SummarizedExperiment:
#' Since [bcbioRNASeq] is an extension of [SummarizedExperiment], this
#' coercion method is very simple. Here we're simply dropping our `@bcbio` slot,
#' which contains raw cellular barcodes and other bcbio-specific metadata.
setAs(
    from = "bcbioRNASeq",
    to = "SummarizedExperiment",
    .coerceToSummarizedExperiment)
