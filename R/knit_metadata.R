#' Knit metadata table
#'
#' @author Michael Steinbaugh
#'
#' @import tibble
#'
#' @param run \code{bcbio-nextgen} run object
#' @param metadata Metadata data frame
#'
#' @export
knit_metadata <- function(run, metadata) {
    metadata[, unique(c("description", run$intgroup))] %>%
        remove_rownames %>%
        kable(caption = "Sample metadata")
}
