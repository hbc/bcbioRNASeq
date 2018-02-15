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
#' colData[["age"]] <- factor(c(14, 30, 14, 30))
#' colData(bcb) <- colData
#' colData(bcb) %>% glimpse()
#'
#' # These internal objects will also get updated
#' bcbio(bcb, "DESeqDataSet") %>% colData() %>% glimpse()
#' assays(bcb)[["rlog"]] %>% colData() %>% glimpse()
#' assays(bcb)[["vst"]] %>% colData() %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname colData
#' @export
setMethod(
    "colData<-",
    signature(
        x = "bcbioRNASeq",
        value = "DataFrame"
    ),
    function(x, ..., value) {
        assert_are_identical(ncol(x), nrow(value))
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
    })
