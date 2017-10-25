#' Check Interesting Groups
#'
#' Prevent unwanted downstream behavior when a missing interesting group
#' is requested by the user.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param object Object supporting [colnames()], typically a [data.frame].
#' @param interestingGroups Interesting groups character vector.
#' @param warnOnNULL Warn the user on `NULL` argument.
#'
#' @return Valid character of defined interesting groups. Stop on failure.
.checkInterestingGroups <- function(
    object,
    interestingGroups,
    warnOnNULL = FALSE) {
    # Check that interesting groups are present in DESeqTransform colData
    if(!all(interestingGroups %in% colnames(object))) {
        stop(paste(
            "Interesting group(s) not defined in metadata:",
            toString(setdiff(interestingGroups, colnames(object)))
        ), call. = FALSE)
    }
    # Default to `sampleName` if `NULL`
    if (is.null(interestingGroups)) {
        if (isTRUE(warnOnNULL)) {
            warning(paste(
                "'interestingGroups' is 'NULL'.",
                "Defaulting to 'sampleName'."
            ), call. = FALSE)
        }
        interestingGroups <- "sampleName"
    }
    interestingGroups
}
