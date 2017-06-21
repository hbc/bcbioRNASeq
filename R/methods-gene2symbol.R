#' Gene identifier to symbol conversion
#'
#' @author Michael Steinbaugh
#'
#' @rdname gene2symbol
#' @docType methods
#'
#' @param object Object.
#'
#' @return [data.frame].
#' @export
setMethod("gene2symbol", "bcbioRnaDataSet", function(object) {
    rowData(bcb)[, c("ensgene", "symbol")] %>%
        as.data.frame %>%
        remove_na %>%
        set_rownames(.[[1L]])
})
