#' Get Flat Files from S4 Object
#'
#' Prepare an unstructured list of information saved in the S4 object,
#' for improved archival data storage.
#'
#' @rdname flatFiles
#' @name flatFiles
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom bcbioBase flatFiles
#'
#' @inheritParams general
#'
#' @return [list].
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' flatFiles(bcb) %>% names()
NULL



# Constructors =================================================================
.flatFiles.bcbioRNASeq <- function(object) {  # nolint
    list(
        assays = assays(object),
        rowData = rowData(object),
        colData = colData(object),
        metadata = metadata(object),
        bcbio = bcbio(object)
    )
}



# Methods ======================================================================
#' @rdname flatFiles
#' @export
setMethod(
    "flatFiles",
    signature("bcbioRNASeq"),
    .flatFiles.bcbioRNASeq)
