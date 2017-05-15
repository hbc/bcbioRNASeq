#' Accessors for the count matrix of a \linkS4class{bcbioRnaDataSet} object
#'
#' @rdname txi
#' @keywords internal
#' @docType methods
#' @name txi
#' @aliases txi txi,bcbioRnaDataSet-method txi<-,bcbioRnaDataSet,matrix-method
#'
#' @param object \linkS4class{bcbioRnaDataSet} object.
#' @param value An integer matrix.
#' @param ... Matrix count data.
#'
#' @return [base::matrix] with raw count data.



#' @rdname txi
#' @export
setGeneric("txi", function(object, ...) standardGeneric("txi"))

#' @rdname txi
#' @export
setGeneric("txi<-", function(object, ..., value) standardGeneric("txi<-"))




#' @rdname txi
#' @export
txi.bcbioRnaDataSet <- function(object) {
    assays(object)[["counts"]]
}



#' @rdname txi
#' @export
setMethod("txi",
          signature(object = "bcbioRnaDataSet"),
          txi.bcbioRnaDataSet)



#' @rdname txi
#' @name "<-txi"
#' @exportMethod "txi<-"
setReplaceMethod(
    "txi", signature(object = "bcbioRnaDataSet",
                     value = "matrix"),
    function(object, value) {
        assays(object)[["counts"]] <- value
        validObject(object)
        object
    })
