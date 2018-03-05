#' Column Data
#'
#' @description
#' Improved assignment method support for [bcbioRNASeq] object.
#'
#' This method support will also update the `colData` inside the `bcbio` and
#' `assays` slots.
#'
#' @rdname colData
#' @name colData
#'
#' @inheritParams general
#'
#' @seealso
#' `help("colData", "SummarizedExperiment")`
#'
#' @return `DataFrame`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # Return as data.frame
#' colData(bcb, return = "data.frame")
#'
#' # Assignment support
#' colData <- colData(bcb)
#' # All columns will be coerced to factors
#' colData[["age"]] <- c(14L, 30L, 14L, 30L)
#' colData(bcb) <- colData
#' colData(bcb) %>% glimpse()
#'
#' # These internal objects will also get updated
#' bcbio(bcb, "DESeqDataSet") %>% colData() %>% glimpse()
#' assays(bcb)[["rlog"]] %>% colData() %>% glimpse()
#' assays(bcb)[["vst"]] %>% colData() %>% glimpse()
NULL



# Constructors =================================================================
.colData <- function(x, return = c("DataFrame", "data.frame")) {
    return <- match.arg(return)
    slot(x, "colData") %>%
        as(return)
}



#' @importFrom basejump sanitizeColData
`.colData<-` <- function(x, ..., value) {
    assert_are_identical(colnames(x), rownames(value))
    value <- as(value, "DataFrame")

    # Sanitize all columns as factors
    value <- sanitizeColData(value)
    assert_has_dimnames(value)

    if (!is.null(bcbio(x, "DESeqDataSet"))) {
        colData(bcbio(x, "DESeqDataSet")) <- value
    }
    if (!is.null(assays(x)[["rlog"]])) {
        colData(assays(x)[["rlog"]]) <- value
    }
    if (!is.null(assays(x)[["vst"]])) {
        colData(assays(x)[["vst"]]) <- value
    }

    slot(x, "colData") <- value
    x
}



# Methods ======================================================================
#' @rdname colData
#' @export
setMethod(
    "colData",
    signature("bcbioRNASeq"),
    .colData
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
