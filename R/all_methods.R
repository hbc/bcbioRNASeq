#' Accessors for the count matrix of a [bcbioRnaDataSet] object
#'
#' This method will be used to access all different count matrix from
#' the object. Gene expression, transcript expression, miRNA expression...
#'
#' @rdname bcbio
#' @docType methods
#' @name bcbio
#' @aliases bcbio bcbio,bcbioRnaDataSet-method bcbio<-,bcbioRnaDataSet, bcbioNames
#'
#' @param object [bcbioRnaDataSet] object.
#' @param value An integer matrix or other object.
#' @param type type of count data to retrieve
#' @param ... Matrix count data.
#'
#' @return Matrix/Object containing count data.

#' @rdname bcbio
#' @export
setGeneric("bcbioNames", function(object, ...) standardGeneric("bcbioNames"))

#' @rdname bcbio
#' @export
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))

#' @rdname bcbio
#' @export
bcbio.bcbioRnaDataSet <- function(object, type="counts") {
    if (type == "counts")
        return(assays(object)[["counts"]])
    if (type %in% names(slot(object, "callers")))
        return(slot(object, "callers")[[type]])
    message(type, " not found.")
}

bcbio.names <- function(object){
    names(slot(object, "callers"))
}

#' @rdname bcbio
#' @export
setMethod("bcbio",
          signature(object = "bcbioRnaDataSet"),
          bcbio.bcbioRnaDataSet)

#' @rdname bcbio
#' @export
setMethod("bcbioNames",
          signature(object = "bcbioRnaDataSet"),
          bcbio.names)

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
            slot(object, "callers")[[type]] <- value
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
