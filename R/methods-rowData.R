#' Row Data
#'
#' @name rowData
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param return Return as "`data.frame`" or "`DataFrame`".
#'
#' @return Data describing the rows of the object.
#'
#' @examples
#' rowData(bcb_small) %>% glimpse()
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
