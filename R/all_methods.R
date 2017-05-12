#' @rdname txi
#' @export
setGeneric("txi", function(object, ...) standardGeneric("txi"))
#' @rdname txi
#' @export
setGeneric("txi<-", function(object, ..., value) standardGeneric("txi<-"))

#' Accessors for the count matrix of a bcbioRnaDataSet object.
#'
#'
#' @docType methods
#' @name txi
#' @rdname txi
#' @aliases txi txi,bcbioRnaDataSet-method txi<-,bcbioRnaDataSet,matrix-method
#'
#' @param object a \code{bcbioRnaDataSet} object
#' @param value an integer matrix
#' @param ... matrix count data
#'
#' @return \code{\link[base]{matrix}} with raw count data.
#' @author Lorena Pantano
#' @export
txi.bcbioRnaDataSet <- function(object) {
    assays(object)[["counts"]]
}

#' @rdname txi
#' @export
setMethod("txi", signature(object="bcbioRnaDataSet"), txi.bcbioRnaDataSet)

#' @name "<-txi"
#' @rdname txi
#' @exportMethod "txi<-"
setReplaceMethod("txi", signature(object="bcbioRnaDataSet", value="matrix"),
                 function(object, value){
                     assays(object)[["counts"]] <- value
                     validObject(object)
                     object
                 }
)
