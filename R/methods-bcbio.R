#' [bcbioRnaDataSet] count matrix accessors
#'
#' This method will be used to access all different count matrix from
#' the object. Gene expression, transcript expression, miRNA expression...
#'
#' @rdname bcbio
#' @docType methods
#'
#' @param object [bcbioRnaDataSet] object.
#' @param type Type of count data to retrieve.
#' @param value An integer matrix or other object.
#' @param ... Additional arguments.
#'
#' @return Count matrix.
#' @export
setMethod("bcbio", "bcbioRnaDataSet", function(object, type = "counts") {
    if (!type %in% count_slots) {
        stop("Unsupported type")
    }
    assays(object)[[type]]
})

#' @rdname bcbio
#' @exportMethod "bcbio<-"
setReplaceMethod(
    "bcbio",
    signature(object = "bcbioRnaDataSet", value = "matrix"),  # ANY
    function(object, type = "counts", value) {
        if (!type %in% c("counts", "abundance", "length", "alt_counts")) {
            stop("Unsupported type")
        }
        assays(object)[[type]] <- value
        validObject(object)
        object
    })
