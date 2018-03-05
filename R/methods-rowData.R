#' Row Data
#'
#' @name rowData
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @inheritParams general
#'
#' @seealso `help("rowData", "SummarizedExperiment")`
#'
#' @return `DataFrame` or `data.frame`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # Return as data.frame
#' rowData(bcb, return = "data.frame")
NULL



# Constructors =================================================================
.rowData <- function(x, return = c("DataFrame", "data.frame", "GRanges")) {
    return <- match.arg(return)
    data <- slot(x, "elementMetadata")
    names <- slot(x, "NAMES")

    data <- as(data, return)

    if (has_dims(data)) {
        rownames(data) <- names
    } else if (has_names(data)) {
        names(data) <- names
    }

    data
}
