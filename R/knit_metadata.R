#' Knit metadata table
#'
#' @author Michael Steinbaugh
#'
#' @importFrom tibble remove_rownames
#'
#' @param run \code{bcbio-nextgen} run
#' @param metadata Metadata data frame
#'
#' @export
knit_metadata <- function(run, metadata) {
    metadata[, unique(c("description", run$intgroup))] %>%
        remove_rownames %>%
        kable(caption = "Sample metadata")
}
