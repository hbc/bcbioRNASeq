.melt_log10 <- function(counts) {
    counts %>%
        as.data.frame %>%
        rownames_to_column %>%
        melt(id = 1L) %>%
        set_names(c("rowname",  # ensembl_gene_id
                    "colname",  # description
                    "counts")) %>%
        # Filter zero counts
        filter(.data[["counts"]] > 0L) %>%
        # log10 transform
        mutate(counts = log10(.data[["counts"]]),
               # [melt()] sets colnames as factor
               colname = as.character(.data[["colname"]]))
}



.join_melt <- function(counts, metadata) {
    if (!identical(colnames(counts), metadata[["colname"]])) {
        stop("Sample description mismatch between counts and metadata")
    }
    .melt_log10(counts) %>%
        left_join(metadata, by = "colname") %>%
        rename(description = .data[["colname"]])
}



#' Melt count matrix to long format and log10 transform
#'
#' @rdname melt_log10
#' @docType methods
#' @author Michael Steinbaugh
#'
#' @param object Object containing counts matrix.
#' @param normalized Select normalized counts (`TRUE`), raw counts (`FALSE`),
#' or specifically request TMM-normalized counts (`tmm`).
#' @param interesting_groups *Optional*. Interesting groups.
#' @param ... Optional parameters.
#'
#' @return log10 melted [data.frame].
#' @export
#'
#' @seealso [reshape2::melt()].
setMethod("melt_log10", "bcbioRNADataSet", function(
    object,
    normalized = TRUE,
    interesting_groups = NULL) {
    counts <- counts(object, normalized = normalized)

    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(object)[["interesting_groups"]]
    }

    metadata <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column("colname") %>%
        tidy_select(!!!syms(c("colname", interesting_groups)))

    .join_melt(counts, metadata)
})

#' @rdname melt_log10
#' @export
setMethod("melt_log10", "DESeqDataSet", function(
    object,
    normalized = FALSE,
    interesting_groups = NULL) {
    counts <- counts(object, normalized = normalized)

    metadata <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column("colname")

    if (!is.null(interesting_groups)) {
        metadata <- metadata %>%
            tidy_select(!!!syms(c("colname", interesting_groups)))
    }

    .join_melt(counts, metadata)
})



#' @rdname melt_log10
#' @export
setMethod("melt_log10", "DESeqTransform", function(
    object,
    interesting_groups = NULL) {
    counts <- assay(object)

    metadata <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column("colname")

    if (!is.null(interesting_groups)) {
        metadata <- metadata %>%
            tidy_select(!!!syms(c("colname", interesting_groups)))
    }

    .join_melt(counts, metadata)
})
