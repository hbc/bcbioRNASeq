#' Unite Interesting Groups
#'
#' Create a single interesting groups column (`interestingGroups`) used for
#' coloring in plots. When multiple interesting groups are present, unite into a
#' single column, delimited by ` : `.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param object Object (e.g. [data.frame]) containing interesting groups
#'   columns.
#' @param interestingGroups Character vector of interesting groups.
#'
#' @return [data.frame].
.uniteInterestingGroups <- function(object, interestingGroups) {
    # Set up the interesting groups column
    interestingGroups <- .checkInterestingGroups(object, interestingGroups)
    if (length(interestingGroups) > 1) {
        object <- unite(
            data = object,
            col = interestingGroups,
            !!!syms(interestingGroups),
            sep = ":",
            remove = FALSE)
    } else {
        object[["interestingGroups"]] <- object[[interestingGroups]]
    }
    object
}
