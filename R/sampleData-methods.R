#' @name sampleData
#' @author Michael Steinbaugh
#' @inherit basejump::sampleData
#' @inheritParams params
#'
#' @param clean `logical(1)`.
#'   Only return `factor` columns not defined in
#'   [metadataBlacklist][bcbioBase::metadataBlacklist]. This removes metrics
#'   columns used for quality control analysis, which are often not informative
#'   as sample metadata.
#'
#' @examples
#' data(bcb)
#' sampleData(bcb, clean = TRUE) %>% lapply(head)
#'
#' ## Assignment support.
#' x <- bcb
#' sampleData(x)[["test"]] <- seq_len(ncol(x))
#' ## `test` column should be now defined.
#' "test" %in% colnames(sampleData(x))
NULL



#' @importFrom basejump sampleData
#' @aliases NULL
#' @export
basejump::sampleData



# FIXME Move `metadataBlacklist` to basejump.
sampleData.bcbioRNASeq <-  # nolint
    function(object, clean = FALSE) {
        assert(isFlag(clean))
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
    definition = sampleData.bcbioRNASeq
)
