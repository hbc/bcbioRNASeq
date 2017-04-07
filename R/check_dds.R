#' @rdname check
#' @param dds DESeqDataSet
#' @param stop Stop if class doesn't match
#' @export
check_dds <- function(dds, stop = TRUE) {
    if (class(dds)[1] == "DESeqDataSet") {
        return(TRUE)
    } else {
        if (isTRUE(stop)) {
            stop("DESeqDataSet required")
        } else {
            return(FALSE)
        }
    }
}
