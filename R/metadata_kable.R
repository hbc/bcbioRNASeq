#' Metadata kable
#'
#' Returns a subset of metadata columns of interest used for knit reports. These
#' "interesting group" columns are defined as \code{intgroup} in the
#' \code{bcbio-nextgen} \code{run} object.
#'
#' @author Michael Steinbaugh
#'
#' @importFrom tibble remove_rownames
#'
#' @param run \code{bcbio-nextgen} run
#' @param metadata Metadata data frame
#'
#' @return Data frame containing only the columns of interest
#' @export
metadata_kable <- function(run, metadata) {
    metadata[, unique(c("description", run$intgroup))] %>%
        remove_rownames %>%
        kable(caption = "Sample metadata")
}
