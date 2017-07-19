#' @rdname aggregate_replicates
#' @usage NULL
#' @export
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})



#' @rdname bcbio
#' @usage NULL
#' @export
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))



#' @rdname bcbio
#' @usage NULL
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))



#' @rdname melt_log10
#' @usage NULL
#' @export
setGeneric("melt_log10", function(object, ...) standardGeneric("melt_log10"))



#' @rdname metadata_table
#' @usage NULL
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})



#' @rdname metrics
#' @usage NULL
#' @export
setGeneric("metrics", function(object) standardGeneric("metrics"))



#' @rdname sample_dirs
#' @usage NULL
#' @export
setGeneric("sample_dirs", function(object) standardGeneric("sample_dirs"))



#' @rdname tmm
#' @usage NULL
#' @export
setGeneric("tmm", function(object) standardGeneric("tmm"))



#' @rdname tpm
#' @usage NULL
#' @export
setGeneric("tpm", function(object) standardGeneric("tpm"))



#' @rdname txi
#' @usage NULL
#' @export
setGeneric("txi", function(object) standardGeneric("txi"))
