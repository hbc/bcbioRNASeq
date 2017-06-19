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
#' Selects `groups_of_interest`-defined interesting groups automatically.
#'
#' @rdname select_groups_of_interest
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param dds [DESeqDataSet] or [DESeqTransform].
#' @param groups_of_interest Groups of interest.
#'
#' @return [colData()] data frame.
.select_groups_of_interest <- function(dds, groups_of_interest) {
    colData(dds) %>%
        as.data.frame %>%
        tidy_select(!!quo(groups_of_interest)) %>%
        DataFrame
}
