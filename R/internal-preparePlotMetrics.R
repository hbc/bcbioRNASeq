#' Prepare Plot Metrics Constructor
#'
#' Define the interesting groups column used for coloring in quality control
#' plots. When multiple interesting groups are present, unite into a single
#' column, delimited by ` : `.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param metrics [metrics()] return.
#' @param interestingGroups [interestingGroups()] return.
#'
#' @return [data.frame].
.preparePlotMetrics <- function(metrics, interestingGroups) {
    # Set up the interesting groups column
    interestingGroups <- .checkInterestingGroups(metrics, interestingGroups)
    if (length(interestingGroups) > 1) {
        metrics <- unite(
            data = metrics,
            col = interestingGroups,
            !!!syms(interestingGroups),
            sep = " : ",
            remove = FALSE)
    } else {
        metrics[["interestingGroups"]] <- metrics[[interestingGroups]]
    }
    metrics
}
