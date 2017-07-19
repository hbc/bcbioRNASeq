#' Metadata Table
#'
#' Returns a subset of metadata columns of interest used for knit reports. These
#' "interesting group" columns are defined as `interesting_groups` in the
#' [bcbioRNADataSet] object.
#'
#' @rdname metadata_table
#' @author Michael Steinbaugh
#'
#' @param ... Additional parameters, passed to [kable()].
#'
#' @return [data.frame] containing only the columns of interest.



#' @rdname metadata_table
#' @usage NULL
.metadata_table <- function(object, ...) {
    object %>%
        colData %>%
        as.data.frame %>%
        remove_rownames %>%
        kable(caption = "Sample metadata", ...)
}



#' @rdname metadata_table
#' @export
setMethod("metadata_table", "bcbioRNADataSet", function(object, ...) {
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
