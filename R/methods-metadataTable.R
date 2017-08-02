#' Metadata Table
#'
#' Returns a subset of metadata columns of interest used for knit reports. These
#' "interesting group" columns are defined as `interestingGroups` in the
#' [bcbioRNADataSet] object.
#'
#' @rdname metadataTable
#' @author Michael Steinbaugh
#'
#' @param ... Additional parameters, passed to [kable()].
#'
#' @return [data.frame] containing only the columns of interest.



#' @rdname metadataTable
#' @usage NULL
.metadataTable <- function(object, ...) {
    object %>%
        colData %>%
        as.data.frame %>%
        remove_rownames %>%
        kable(caption = "Sample metadata", ...)
}



#' @rdname metadataTable
#' @export
setMethod("metadataTable", "bcbioRNADataSet", .metadataTable)



#' @rdname metadataTable
#' @export
setMethod("metadataTable", "DESeqDataSet", .metadataTable)



#' @rdname metadataTable
#' @export
setMethod("metadataTable", "DESeqTransform", .metadataTable)
