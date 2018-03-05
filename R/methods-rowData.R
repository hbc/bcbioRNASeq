#' Row Data
#'
#' @name rowData
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @inheritParams general
#'
#' @param return Return as "`DataFrame`", "`data.frame`", or unmodified
#'   ("`AsIs`").
#'
#' @return Data describing the rows of the object.
#'
#' @seealso `help("rowData", "SummarizedExperiment")`
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # Return as data.frame
#' rowData(bcb, return = "data.frame")
NULL



# Constructors =================================================================
.rowData <- function(x, return = c("DataFrame", "data.frame", "AsIs")) {
    return <- match.arg(return)
    data <- slot(x, "elementMetadata")
    if (return != "AsIs") {
        data <- as(data, return)
    }
    names <- slot(x, "NAMES")
    if (has_dims(data)) {
        rownames(data) <- names
    } else if (has_names(data)) {
        names(data) <- names
    }
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
