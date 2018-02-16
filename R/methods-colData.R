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
#' @return [DataFrame].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
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
`.colData<-` <- function(x, ..., value) {
    assert_are_identical(colnames(x), rownames(value))

    # Sanitize all columns as factors
    value <- lapply(
        X = value,
        FUN = function(x) {
            droplevels(as.factor(x))
        }
    ) %>%
        as("DataFrame")

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



# Assignment Methods ===========================================================
#' @rdname colData
#' @export
setMethod(
    "colData<-",
    signature(
        x = "bcbioRNASeq",
        value = "data.frame"
    ),
    `.colData<-`)



#' @rdname colData
#' @export
setMethod(
    "colData<-",
    signature(
        x = "bcbioRNASeq",
        value = "DataFrame"
    ),
    `.colData<-`)
