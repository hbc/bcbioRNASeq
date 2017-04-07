#' Knit metadata table
#'
#' @author Michael Steinbaugh
#'
#' @param bcbio \code{bcbio-nextgen} run object
#' @param metadata Metadata data frame
#'
#' @export
knit_metadata <- function(bcbio, metadata) {
    metadata[, unique(c("description", bcbio$intgroup))] %>%
        kable(caption = "Sample metadata")
}
