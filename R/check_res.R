check_res <- function(res, stop = TRUE) {
    if (class(res)[1] == "DESeqResults") {
        return(TRUE)
    } else {
        if (isTRUE(stop)) {
            stop("DESeqDataSet required")
        } else {
            return(FALSE)
        }
    }
}
