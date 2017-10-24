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
#'
#' @return Stop on failure.
.checkInterestingGroups <- function(object, interestingGroups) {
    # Check that interesting groups are present in DESeqTransform colData
    if(!all(interestingGroups %in% colnames(object))) {
        stop(paste(
            "Missing interesting groups:",
            toString(setdiff(interestingGroups, colnames(object)))
        ), call. = FALSE)
    }
}
