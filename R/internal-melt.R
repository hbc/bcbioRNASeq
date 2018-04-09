#' Melt Counts Matrix to Long Format
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @inheritParams general
#' @param normalized Select normalized counts (`TRUE`), raw counts (`FALSE`),
#' or specifically request TMM-normalized counts (`tmm`).
#'
#' @seealso [reshape2::melt()].
#'
#' @return Melted `grouped_df`, grouped by `sampleID`.
#'
#' @examples
#' counts <- counts(bcb_small)
#' sampleData <- sampleData(bcb_small)
#' .meltCounts(counts, sampleData)
.meltCounts <- function(object, sampleData = NULL) {
    assert_is_matrix(object)

    data <- object %>%
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
            colnames(object),
            as.character(sampleData[["sampleID"]])
        )
        data <- left_join(
            x = data,
            y = as.data.frame(sampleData),
            by = "sampleID"
        )
    }

    data
}



.meltLog2Counts <- function(object, ...) {
    x <- .meltCounts(object, ...)
    x[["counts"]] <- log2(x[["counts"]] + 1L)
    x
}
