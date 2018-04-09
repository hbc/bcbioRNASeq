#' Melt Counts Matrix to Long Format
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @seealso [reshape2::melt()].
#'
#' @return Melted `grouped_df`, grouped by `sampleID` and `geneID`.
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

    if (!is.null(sampleData)) {
        assert_has_dims(sampleData)
        assert_are_identical(
            colnames(counts),
            as.character(sampleData[["sampleID"]])
        )
        data <- left_join(
            x = data,
            y = as.data.frame(sampleData),
            by = "sampleID"
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
