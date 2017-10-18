#' Get Flat Files from S4 Object
#'
#' Prepare an unstructured list of information saved in the S4 object,
#' for improved archival data storage.
#'
#' @rdname flatFiles
#' @name flatFiles
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams AllGenerics
#'
#' @return [list].
NULL



# Methods ====
#' @rdname flatFiles
#' @export
setMethod(
    "flatFiles",
    signature("bcbioRNASeqANY"),
    function(object) {
        list(
            assays = assays(object),
            rowData = rowData(object),
            colData = colData(object),
            metadata = metadata(object),
            bcbio = bcbio(object)
        )
    })
