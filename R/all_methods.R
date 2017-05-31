#' Accessors for the count matrix of a [bcbioRnaDataSet] object
#'
#' @rdname txi
#' @keywords internal
#' @docType methods
#' @name txi
#' @aliases txi txi,bcbioRnaDataSet-method txi<-,bcbioRnaDataSet,matrix-method
#'
#' @param object [bcbioRnaDataSet] object.
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


#' Accessors for the count matrix of a [bcbioRnaDataSet] object
#'
#' This method will be used to access all different count matrix from
#' the object. Gene expression, transcript expression, miRNA expression...
#'
#' @rdname bcbio
#' @docType methods
#' @name bcbio
#' @aliases bcbio bcbio,bcbioRnaDataSet-method bcbio<-,bcbioRnaDataSet,matrix-method
#'
#' @param object [bcbioRnaDataSet] object.
#' @param value An integer matrix or other object.
#' @param type type of count data to retrieve
#' @param ... Matrix count data.
#'
#' @return Matrix/Object containing count data.

#' @rdname bcbio
#' @export
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))

#' @rdname bcbio
#' @export
bcbio.bcbioRnaDataSet <- function(object, type="counts") {
    if (type == "counts")
        return(assays(object)[["counts"]])
    if (type %in% names(metadata(object)[["alt_counts"]]))
        return(metadata(object)[["alt_counts"]][[type]])
    message(type, " not found.")
}

#' @rdname bcbio
#' @export
setMethod("bcbio",
          signature(object = "bcbioRnaDataSet"),
          bcbio.bcbioRnaDataSet)

#' @rdname bcbio
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))

#' @rdname bcbio
#' @exportMethod "bcbio<-"
setReplaceMethod(
    "bcbio", signature(object = "bcbioRnaDataSet",
                      value = "ANY"),
    function(object, type="counts", value) {
        if (type == "counts"){
            assays(object)[["counts"]] <- value
            validObject(object)
        }else{
            metadata(object)[["alt_counts"]][[type]] <- value
        }
        object
    })


bcbio.samples <- function(object){
    metadata(object)[["sample_dirs"]]
}
#' Accessors samples dir of a [bcbioRnaDataSet] object
#'
#' This method will be used to access folders where
#' sample information is kept
#'
#' @rdname bcbcols
#' @docType methods
#' @name bcbcols
#' @aliases bcbcols bcbcols,bcbioRnaDataSet
#'
#' @param object [bcbioRnaDataSet] object.
#'
#' @return folders where samples are kept
setGeneric("bcbcols", function(object) standardGeneric("bcbcols"))
#' @rdname bcbcols
#' @export
setMethod("bcbcols",
          signature(object = "bcbioRnaDataSet"),
          bcbio.samples)
