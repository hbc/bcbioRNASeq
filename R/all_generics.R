#' aggregate_replicates ====
#' Aggregate lane-split technical replicates
#'
#' Frequently RNA-seq experiments are performed with technical replicates
#' split across flow cell lanes. This generic facilitates quick aggregation
#' of counts across the flow cells.
#'
#' @rdname aggregate_replicates
#' @name aggregate_replicates
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param ... Passthrough.
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})



# bcbio ====
#' [bcbioRnaDataSet] count matrix accessors
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
#' @return Count matrix.
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))

#' @rdname bcbio
#' @name "bcbio<-"
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))



# melt_log10 ====
#' Melt count matrix to long format and log10 transform
#'
#' @rdname melt_log10
#' @name melt_log10
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Object containing counts matrix.
#' @param ... Optional parameters.
#'
#' @return log10 melted data frame.
#' @export
setGeneric("melt_log10", function(object, ...) standardGeneric("melt_log10"))



# metadata_table ====
#' Metadata table
#'
#' Returns a subset of metadata columns of interest used for knit reports. These
#' "interesting group" columns are defined as `interesting_groups` in the
#' [bcbioRnaDataSet] object.
#'
#' @rdname metadata_table
#' @name metadata_table
#' @docType methods
#'
#' @param ... [kable()] passthrough parameters.
#'
#' @return Data frame containing only the columns of interest.
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})


# metrics ====
#' Sample metrics
#'
#' @rdname metrics
#' @name metrics
#' @docType methods
#'
#' @param object [bcbioRnaDataSet] object.
#'
#' @export
setGeneric("metrics", function(object) standardGeneric("metrics"))



# sample_dirs ====
#' Sample directories of [bcbioRnaDataSet] object.
#'
#' This method will be used to access folders where sample information is kept.
#'
#' @rdname sample_dirs
#' @name sample_dirs
#' @docType methods
#'
#' @param object [bcbioRnaDataSet] object.
#'
#' @return Folders where samples are kept.
#' @export
setGeneric("sample_dirs", function(object) standardGeneric("sample_dirs"))



# tmm ===
#' Trimmed mean of M-values (TMM) normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @author Michael Steinbaugh
#'
#' @param object Object containing counts.
#'
#' @export
setGeneric("tmm", function(object) standardGeneric("tmm"))



# tpm ====
#' Transcripts per million (TPM)
#'
#' @rdname tpm
#' @name tpm
#' @docType methods
#'
#' @param object [bcbioRnaDataSet] object.
#'
#' @export
setGeneric("tpm", function(object) standardGeneric("tpm"))
