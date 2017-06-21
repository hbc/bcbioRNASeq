#' @rdname aggregate_replicates
#' @export
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})

#' @rdname bcbio
#' @export
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))

#' @rdname bcbio
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))

#' @rdname melt_log10
#' @export
setGeneric("melt_log10", function(object, ...) standardGeneric("melt_log10"))

#' @rdname metadata_table
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})

#' @rdname metrics
#' @export
setGeneric("metrics", function(object) standardGeneric("metrics"))

#' @rdname sample_dirs
#' @export
setGeneric("sample_dirs", function(object) standardGeneric("sample_dirs"))

#' @rdname tmm
#' @export
setGeneric("tmm", function(object) standardGeneric("tmm"))

#' @rdname tpm
#' @export
setGeneric("tpm", function(object) standardGeneric("tpm"))
