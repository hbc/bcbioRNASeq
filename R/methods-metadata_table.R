.metadata_table <- function(object, ...) {
    object %>%
        colData %>%
        as.data.frame %>%
        remove_rownames %>%
        kable(caption = "Sample metadata", ...)
}

#' Metadata table
#'
#' Returns a subset of metadata columns of interest used for knit reports. These
#' "interesting group" columns are defined as `interesting_groups` in the
#' [bcbioRnaDataSet] object.
#'
#' @rdname metadata_table
#' @docType methods
#'
#' @param object Object.
#' @param ... Additional parameters passed to [kable()].
#'
#' @return Data frame containing only the columns of interest.
#' @export
setMethod("metadata_table", "bcbioRnaDataSet", function(object, ...) {
    .metadata_table(object)
})

#' @rdname metadata_table
#' @export
setMethod("metadata_table", "DESeqDataSet", function(object, ...) {
    .metadata_table(object)
})

#' @rdname metadata_table
#' @export
setMethod("metadata_table", "DESeqTransform", function(object, ...) {
    .metadata_table(object)
})
