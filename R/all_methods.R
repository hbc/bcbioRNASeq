# Count matrix accessors =======================================================
#' Accessors for the count matrix of a [bcbioRnaDataSet] object
#'
#' This method will be used to access all different count matrix from
#' the object. Gene expression, transcript expression, miRNA expression...
#'
#' @rdname bcbio
#' @docType methods
#' @name bcbio
#' @aliases bcbio bcbio,bcbioRnaDataSet-method
#'   bcbio<-,bcbioRnaDataSet,matrix-method
#'
#' @param object [bcbioRnaDataSet] object.
#' @param value An integer matrix or other object.
#' @param type type of count data to retrieve
#' @param ... Matrix count data.
#'
#' @return Matrix/Object containing count data.
NULL



#' @rdname bcbio
#' @export
bcbio.bcbioRnaDataSet <- function(object, type="counts") {  # nolint
    if (type == "counts")
        return(assays(object)[["counts"]])
    if (type %in% names(slot(object, "callers")))
        return(slot(object, "callers")[[type]])
    message(type, " not found.")
}

#' @rdname bcbio
#' @export
setMethod("bcbio",
          signature(object = "bcbioRnaDataSet"),
          bcbio.bcbioRnaDataSet)

#' @rdname bcbio
#' @exportMethod "bcbio<-"
setReplaceMethod(
    "bcbio", signature(object = "bcbioRnaDataSet", value = "ANY"),
    function(object, type="counts", value) {
        if (type == "counts") {
            assays(object)[["counts"]] <- value
            validObject(object)
        }else{
            slot(object, "callers")[[type]] <- value
        }
        object
    })






# Sample accessors =============================================================
bcbio.samples <- function(object){  # nolint
    metadata(object)[["sample_dirs"]]
}



#' Accessors samples dir of a [bcbioRnaDataSet] object.
#'
#' This method will be used to access folders where sample information is kept.
#'
#' @rdname bcbcols
#' @docType methods
#' @name bcbcols
#' @aliases bcbcols bcbcols,bcbioRnaDataSet
#'
#' @param object [bcbioRnaDataSet] object.
#'
#' @return Folders where samples are kept.
#' @export
setMethod("bcbcols", signature(object = "bcbioRnaDataSet"), bcbio.samples)
