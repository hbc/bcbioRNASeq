#' Column Data
#'
#' Improved assignment method support for a `bcbioRNASeq` object.
#'
#' This method support will also update the nested [colData()] inside the
#' `assays` slot.
#'
#' @name colData
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @inheritParams general
#' @param return Return as "`data.frame`" or "`DataFrame`".
#'
#' @return Data describing the columns of the object.
#'
#' @seealso `help("colData", "SummarizedExperiment")`
#'
#' @examples
#' colData(bcb_small) %>% glimpse()
#'
#' # Assignment method support
#' # All columns will be coerced to factors
#' colData(bcb_small)[["age"]] <- c(14L, 30L, 14L, 30L, 14L, 30L)
#' colData(bcb_small) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom basejump sanitizeSampleData
`.colData<-` <- function(x, value) {
    validObject(x)
    assert_are_identical(colnames(x), rownames(value))
    value <- as(value, "DataFrame")
    value <- sanitizeSampleData(value)
    assert_has_dimnames(value)
    slot(x, "colData") <- value
    x
}



# Methods ======================================================================
#' @rdname colData
#' @importFrom basejump sanitizeSampleData
#' @export
setMethod(
    "colData",
    signature("bcbioRNASeq"),
    function(x, return = c("data.frame", "DataFrame")) {
        return <- match.arg(return)
        data <- slot(x, "colData")
        data <- sanitizeSampleData(data)
        as(data, return)
    }
)



# Assignment Methods ===========================================================
#' @rdname colData
#' @export
setMethod(
    "colData<-",
    signature(
        x = "bcbioRNASeq",
        value = "data.frame"
    ),
    `.colData<-`
)



#' @rdname colData
#' @export
setMethod(
    "colData<-",
    signature(
        x = "bcbioRNASeq",
        value = "DataFrame"
    ),
    `.colData<-`
)
