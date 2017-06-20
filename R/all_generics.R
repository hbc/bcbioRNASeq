# bcbio ====
# [TODO] Too vague
#' Accessors for the count matrix of a [bcbioRnaDataSet] object
#'
#' This method will be used to access all different count matrix from
#' the object. Gene expression, transcript expression, miRNA expression...
#'
#' @rdname bcbio
#' @name bcbio
#' @docType methods
#'
#' @aliases bcbio bcbio,bcbioRnaDataSet-method
#'   bcbio<-,bcbioRnaDataSet,matrix-method
#'
#' @param object [bcbioRnaDataSet] object.
#' @param value An integer matrix or other object.
#' @param type type of count data to retrieve
#' @param ... Matrix count data.
#'
#' @return Matrix/Object containing count data.
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))

#' @rdname bcbio
#' @name "bcbio<-"
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))



# metrics ====
#' Sample metrics
#'
#' @rdname metrics
#' @name metrics
#' @export
setGeneric("metrics", function(object) standardGeneric("metrics"))



# raw_counts ====
#' Raw counts
#'
#' @rdname raw_counts
#' @name raw_counts
#' @export
setGeneric("raw_counts", function(object) standardGeneric("raw_counts"))



# sample_dirs ====
#' Sample directories of [bcbioRnaDataSet] object.
#'
#' This method will be used to access folders where sample information is kept.
#'
#' @rdname sample_dirs
#' @docType methods
#' @name bcbcols
#' @aliases bcbcols bcbcols,bcbioRnaDataSet
#'
#' @param object [bcbioRnaDataSet] object.
#'
#' @return Folders where samples are kept.
#' @export
#' @rdname sample_dirs
#' @name sample_dirs
#' @export
setGeneric("sample_dirs", function(object) standardGeneric("sample_dirs"))



# tpm ====
#' Transcripts per million (TPM)
#'
#' @rdname tpm
#' @name tpm
#' @export
setGeneric("tpm", function(object) standardGeneric("tpm"))
