#' Gene Identifier to Symbol Conversion
#'
#' @rdname gene2symbol
#' @author Michael Steinbaugh
#'
#' @return [data.frame].
#' @export
setMethod("gene2symbol", "bcbioRNADataSet", function(object) {
    rowData(object)[, c("ensgene", "symbol")] %>%
        as.data.frame %>%
        remove_na %>%
        set_rownames(.[[1L]])
})
