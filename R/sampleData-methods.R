#' Sample Data
#'
#' Return the sample metadata. Columns are always sanitized to factor.
#'
#' @note This is a complement to the standard [colData()] function, but improves
#'   support for accessing sample metadata for datasets where multiple items in
#'   the columns map to a single sample (e.g. cells for a single-cell RNA-seq
#'   experiment).
#'
#' @name sampleData
#' @family Data Functions
#' @author Michael Steinbaugh
#' @importFrom basejump sampleData
#' @export
#'
#' @inheritParams general
#' @param clean `boolean`. Only return `factor` columns not defined in
#'   [bcbioBase::metadataBlacklist].
#'
#' @return `DataFrame`.
#'
#' @seealso [bcbioBase::metadataBlacklist].
#'
#' @examples
#' # SummarizedExperiment ====
#' sampleData(bcb_small, clean = TRUE) %>% glimpse()
#' sampleData(bcb_small, clean = FALSE) %>% glimpse()
#'
#' # Assignment support
#' x <- bcb_small
#' sampleData(x)[["test"]] <- seq_len(ncol(x))
#' # `test` column should be now defined
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
