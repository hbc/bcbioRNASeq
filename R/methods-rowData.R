#' Row Data
#'
#' @name rowData
#' @author Michael Steinbaugh
#'
#' @param return Return as "`data.frame`" or "`DataFrame`".
#'
#' @return Data describing the rows of the object.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' rowData(bcb) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname rowData
#' @export
setMethod(
    "rowData",
    signature("bcbioRNASeq"),
    function(x, return = c("data.frame", "DataFrame")) {
        return <- match.arg(return)
        data <- mcols(rowRanges(x))
        rownames(data) <- names(rowRanges(x))
        as(data, return)
    }
)
