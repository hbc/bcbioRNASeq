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
#'
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



#' @rdname sampleData
#' @name sampleData<-
#' @importFrom basejump sampleData<-
#' @export
NULL



#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        clean = FALSE
    ) {
        data <- colData(object)
        assert_is_a_bool(clean)

        # Only return factor columns, if desired
        if (isTRUE(clean)) {
            data <- data[, vapply(data, is.factor, logical(1L)), drop = FALSE]
            # Drop remaining blacklisted columns
            setdiff <- setdiff(colnames(data), bcbioBase::metadataBlacklist)
            data <- data[, setdiff, drop = FALSE]
        } else {
            # Include `interestingGroups` column, if not NULL
            if (missing(interestingGroups)) {
                interestingGroups <- bcbioBase::interestingGroups(object)
            }
            if (length(interestingGroups)) {
                data <- uniteInterestingGroups(data, interestingGroups)
            }
        }

        as(data, "DataFrame")
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "bcbioRNASeq",
        value = "DataFrame"
    ),
    function(object, value) {
        value[["interestingGroups"]] <- NULL
        colData(object) <- value
        object
    }
)
