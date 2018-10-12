#' @name sampleData
#' @author Michael Steinbaugh
#' @importFrom basejump sampleData
#' @inherit basejump::sampleData
#' @export
#'
#' @inheritParams general
#' @param clean `boolean`. Only return `factor` columns not defined in
#'   [bcbioBase::metadataBlacklist]. This removes metrics columns used for
#'   quality control analysis, which are often not informative as sample
#'   metadata.
#'
#' @return `DataFrame`.
#'
#' @examples
#' data(bcb_small)
#'
#' sampleData(bcb_small, clean = TRUE) %>% glimpse()
#' sampleData(bcb_small, clean = FALSE) %>% glimpse()
#'
#' ## Assignment support.
#' x <- bcb_small
#' sampleData(x)[["test"]] <- seq_len(ncol(x))
#' ## `test` column should be now defined
#' "test" %in% colnames(sampleData(x))
NULL



.sampleData.bcbioRNASeq <-  # nolint
    function(object, clean = FALSE) {
        assert_is_a_bool(clean)
        rse <- as(object, "RangedSummarizedExperiment")
        data <- sampleData(rse)
        # Only return whitelisted factor columns, if desired.
        if (isTRUE(clean)) {
            data <- data[, vapply(data, is.factor, logical(1L)), drop = FALSE]
            # Drop remaining blacklisted columns.
            setdiff <- setdiff(colnames(data), metadataBlacklist)
            data <- data[, setdiff, drop = FALSE]
        }
        data
    }



#' @rdname sampleData
#' @export
setMethod(
    f = "sampleData",
    signature = signature("bcbioRNASeq"),
    definition = .sampleData.bcbioRNASeq
)
