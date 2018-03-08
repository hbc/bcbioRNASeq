#' Row Data
#'
#' @name rowData
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @inheritParams general
#'
#' @param return Return as "`DataFrame`" or "`data.frame`".
#'
#' @return Data describing the rows of the object.
#'
#' @seealso `help("rowData", "SummarizedExperiment")`
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' rowData(bcb, return = "DataFrame") %>% glimpse()
#' rowData(bcb, return = "data.frame") %>% glimpse()
NULL



# Constructors =================================================================
.rowData <- function(x, return = c("DataFrame", "data.frame")) {
    return <- match.arg(return)
    data <- slot(x, "elementMetadata")
    data <- as(data, return)
    rownames(data) <- slot(x, "NAMES")
    data
}



# Methods ======================================================================
#' @rdname rowData
#' @export
setMethod(
    "rowData",
    signature("bcbioRNASeq"),
    .rowData
)
