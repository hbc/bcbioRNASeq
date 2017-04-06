check_dt <- function(dt, stop = TRUE) {
    if (class(dt)[1] == "DESeqTransform") {
        return(TRUE)
    } else {
        if (isTRUE(stop)) {
            stop("DESeqDataSet required")
        } else {
            return(FALSE)
        }
    }
}
