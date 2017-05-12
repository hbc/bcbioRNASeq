#' Accessors for the count matrix of a bcbioRnaDataSet object.
#'
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,bcbioRnaDataSet-method counts<-,bcbioRnaDataSet,matrix-method
#'
#' @param object a \code{bcbioRnaDataSet} object
#' @param value an integer matrix
#' @return \code{\link[base]{matrix}} with raw count data.
#' @author Lorena Pantano
#' @export
counts.bcbioRnaDataSet <- function(object) {
    assays(object)[["counts"]]
}

#' @rdname counts
#' @export
setMethod("counts", signature(object="bcbioRnaDataSet"), counts.bcbioRnaDataSet)

#' @name counts
#' @rdname counts
#' @exportMethod "counts<-"
setReplaceMethod("counts", signature(object="bcbioRnaDataSet", value="matrix"),
                 function(object, value){
                     assays(object)[["counts"]] <- value
                     validObject(object)
                     object
                 }
)
