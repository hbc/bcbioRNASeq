#' Melt Count Matrix to Long Format and log10 Transform
#'
#' @rdname meltLog10
#' @name meltLog10
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param normalized Select normalized counts (`TRUE`), raw counts (`FALSE`),
#' or specifically request TMM-normalized counts (`tmm`).
#'
#' @seealso [reshape2::melt()].
#'
#' @return log10 melted [data.frame].
#'
#' @examples
#' data(bcb, dds, rld)
#'
#' # bcbioRNASeq
#' meltLog10(bcb) %>%
#'     str()
#'
#' # DESeqDataSet
#' meltLog10(dds) %>%
#'     str()
#'
#' # DESeqTransform
#' meltLog10(rld) %>%
#'     str()
NULL



# Constructors ====
.joinMelt <- function(counts, metadata) {
    if (!identical(colnames(counts), metadata[["sampleID"]])) {
        stop("Sample name mismatch between counts and metadata")
    }
    .meltLog10(counts) %>%
        left_join(as.data.frame(metadata), by = "sampleID")
}



.meltLog10 <- function(counts) {
    counts %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        melt(id = 1) %>%
        setNames(c("ensgene", "sampleID", "counts")) %>%
        .[.[["counts"]] > 0, ] %>%
        # log10 transform the counts
        mutate(counts = log10(.data[["counts"]]),
               # `melt()` sets colnames as factor
               sampleID = as.character(.data[["sampleID"]]))
}



# Methods ====
#' @rdname meltLog10
#' @export
setMethod("meltLog10", "bcbioRNASeqANY", function(
    object,
    normalized = TRUE) {
    .joinMelt(
        counts = counts(object, normalized = normalized),
        metadata = colData(object))
})



#' @rdname meltLog10
#' @export
setMethod("meltLog10", "DESeqDataSet", function(
    object,
    normalized = TRUE) {
    .joinMelt(
        counts = counts(object, normalized = normalized),
        metadata = colData(object))
})



#' @rdname meltLog10
#' @export
setMethod("meltLog10", "DESeqTransform", function(object) {
    .joinMelt(
        counts = assay(object),
        metadata = colData(object))
})
