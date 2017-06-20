# bcbio ====
#' @rdname bcbio
#' @export
setMethod("bcbio", "bcbioRnaDataSet", function(object, type = "counts") {
    if (!type %in% c("counts", "abundance", "length", "alt_counts")) {
        stop("Unsupported type")
    }
    assays(object)[[type]]
})

# [TODO] Modified to just work on assays (simpler). Otherwise we can define
# our own slots, which may not be necessary.
#' @rdname bcbio
#' @exportMethod "bcbio<-"
setReplaceMethod(
    "bcbio",
    signature(object = "bcbioRnaDataSet", value = "matrix"),  # ANY
    function(object, type = "counts", value) {
        if (!type %in% c("counts", "abundance", "length", "alt_counts")) {
            stop("Unsupported type")
        }
        assays(object)[[type]] <- value
        validObject(object)
        object
    })



# raw_counts ====
#' @rdname raw_counts
#' @export
setMethod("raw_counts", "bcbioRnaDataSet", function(object) {
    assays(object)[["counts"]]
})



# sample_dirs ====
#' @rdname sample_dirs
#' @export
setMethod("sample_dirs", "bcbioRnaDataSet", function(object) {
    metadata(object)[["sample_dirs"]]
})



# tpm ====
#' @rdname tpm
#' @export
setMethod("tpm", "bcbioRnaDataSet", function(object) {
    assays(object)[["abundance"]]
})
