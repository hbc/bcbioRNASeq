#' Melt Counts Matrix to Long Format
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr group_by left_join mutate
#' @importFrom magrittr set_colnames
#' @importFrom reshape2 melt
#' @importFrom rlang !! !!! sym syms
#' @importFrom tibble as_tibble rownames_to_column
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
#' colData <- colData(bcb_small)
#' .meltCounts(counts, colData)
.meltCounts <- function(object, colData) {
    assert_is_matrix(object)
    assert_has_dims(colData)
    assert_are_identical(
        colnames(object),
        as.character(colData[["sampleID"]])
    )
    object %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        melt(id = 1L) %>%
        as_tibble() %>%
        set_colnames(c("geneID", "sampleID", "counts")) %>%
        .[, c("sampleID", "geneID", "counts")] %>%
        arrange(!!!syms(c("sampleID", "geneID"))) %>%
        group_by(!!sym("sampleID")) %>%
        left_join(as.data.frame(colData), by = "sampleID")
}



.meltLog2Counts <- function(object, colData) {
    x <- .meltCounts(object, colData)
    x[["counts"]] <- log2(x[["counts"]] + 1L)
    x
}
