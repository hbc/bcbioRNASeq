#' Melt count matrix to long format
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @seealso [reshape2::melt()].
#'
#' @return `grouped_df`, grouped by `sampleID` and `geneID`.
#'
#' @examples
#' counts <- counts(bcb_small)
#' sampleData <- sampleData(bcb_small)
#' x <- .meltCounts(counts, sampleData)
.meltCounts <- function(counts, sampleData = NULL) {
    assert_is_matrix(counts)

    data <- counts %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        melt(id = 1L) %>%
        as_tibble() %>%
        set_colnames(c("geneID", "sampleID", "counts")) %>%
        arrange(!!!syms(c("sampleID", "geneID"))) %>%
        group_by(!!!syms(c("sampleID", "geneID")))

    if (length(sampleData)) {
        assert_are_set_equal(colnames(counts), rownames(sampleData))
        sampleData[["sampleID"]] <- rownames(sampleData)
        data <- merge(
            x = data,
            y = as.data.frame(sampleData),
            by = "sampleID",
            all.x = TRUE
        )
    }

    if (!"interestingGroups" %in% colnames(data)) {
        data[["interestingGroups"]] <- data[["sampleID"]]
    }

    data
}



.meltLog2Counts <- function(counts, ...) {
    counts <- log2(counts + 1L)
    .meltCounts(counts, ...)
}
