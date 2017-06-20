#' Get contrast name from [DESeqResults]
#'
#' @rdname res_contrast_name
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param res [DESeqResults].
#'
#' @return Contrast name string.
.res_contrast_name <- function(res) {
    mcols(res)[2, 2] %>%
        str_replace("^.*:\\s", "")
}



#' [DESeq] [colData]
#'
#' Selects `interesting_groups`-defined interesting groups automatically.
#'
#' @rdname select_interesting_groups
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param dds [DESeqDataSet] or [DESeqTransform].
#' @param interesting_groups Groups of interest.
#'
#' @return [colData()] data frame.
.select_interesting_groups <- function(dds, interesting_groups) {
    colData(dds) %>%
        as.data.frame %>%
        tidy_select(!!quo(interesting_groups)) %>%
        DataFrame
}
