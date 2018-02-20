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



# Methods ======================================================================
#' @rdname flatFiles
#' @export
setMethod(
    "flatFiles",
    signature("bcbioRNASeq"),
    function(object) {
        list(
            assays = assays(object),
            rowData = rowData(object),
            colData = colData(object),
            metadata = metadata(object),
            bcbio = bcbio(object)
        )
    })
