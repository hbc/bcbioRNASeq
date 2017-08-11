#' Melt Count Matrix to Long Format and log10 Transform
#'
#' @rdname meltLog10
#' @name meltLog10
#'
#' @param normalized Select normalized counts (`TRUE`), raw counts (`FALSE`),
#' or specifically request TMM-normalized counts (`tmm`).
#' @param interestingGroups *Optional*. Interesting groups.
#'
#' @seealso [reshape2::melt()].
#'
#' @return log10 melted [data.frame].
#'
#' @examples
#' data(bcb, dds, rld)
#'
#' # bcbioRNADataSet
#' meltLog10(bcb) %>% glimpse
#'
#' # DESeqDataSet
#' meltLog10(dds) %>% glimpse
#'
#' # DESeqTransform
#' meltLog10(rld) %>% glimpse
NULL



# Methods ====
#' @rdname meltLog10
#' @export
setMethod("meltLog10", "bcbioRNADataSet", function(
    object,
    normalized = TRUE) {
    counts <- counts(object, normalized = normalized)
    interestingGroups <- metadata(object)[["interestingGroups"]]
    metadata <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column("colname") %>%
        .[, c("colname", interestingGroups)]
    .joinMelt(counts, metadata)
})



#' @rdname meltLog10
#' @export
setMethod("meltLog10", "DESeqDataSet", function(
    object,
    normalized = TRUE,
    interestingGroups = NULL) {
    counts <- counts(object, normalized = normalized)
    metadata <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column("colname")
    if (!is.null(interestingGroups)) {
        metadata <- metadata %>%
            .[, c("colname", interestingGroups)]
    }
    .joinMelt(counts, metadata)
})



#' @rdname meltLog10
#' @export
setMethod("meltLog10", "DESeqTransform", function(
    object,
    interestingGroups = NULL) {
    counts <- assay(object)
    metadata <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column("colname")
    if (!is.null(interestingGroups)) {
        metadata <- metadata %>%
            .[, c("colname", interestingGroups)]
    }
    .joinMelt(counts, metadata)
})



#' @rdname meltLog10
#' @usage NULL
.joinMelt <- function(counts, metadata) {
    if (!identical(colnames(counts), metadata[["colname"]])) {
        stop("Sample name mismatch between counts and metadata")
    }
    .meltLog10(counts) %>%
        left_join(metadata, by = "colname") %>%
        rename(sampleName = .data[["colname"]])
}



#' @rdname meltLog10
#' @usage NULL
.meltLog10 <- function(counts) {
    counts %>%
        as.data.frame %>%
        rownames_to_column %>%
        melt(id = 1L) %>%
        setNames(c("rowname",  # ensembl gene ID
                   "colname",  # sample name
                   "counts")) %>%
        .[.[["counts"]] > 0L, ] %>%
        mutate(counts = log10(.data[["counts"]]),
               # [melt()] sets colnames as factor
               colname = as.character(.data[["colname"]]))
}
