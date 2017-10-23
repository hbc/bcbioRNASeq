#' Get Flat Files from S4 Object
#'
#' Prepare an unstructured list of information saved in the S4 object,
#' for improved archival data storage.
#'
#' @rdname flatFiles
#' @name flatFiles
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump flatFiles
#'
#' @inheritParams AllGenerics
#'
#' @return [list].
#'
#' @examples
#' lst <- flatFiles(bcb)
#' names(lst)
NULL



# Methods ====
#' @rdname flatFiles
#' @importFrom S4Vectors metadata
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
