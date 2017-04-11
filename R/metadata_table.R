#' Metadata table
#'
#' Returns a subset of metadata columns of interest used for knit reports. These
#' "interesting group" columns are defined as \code{intgroup} in the
#' \code{bcbio-nextgen} run object.
#'
#' @author Michael Steinbaugh
#'
#' @importFrom knitr kable
#' @importFrom tibble remove_rownames
#'
#' @param run \code{bcbio-nextgen} run
#' @param metadata Metadata data frame
#'
#' @return Data frame containing only the columns of interest
#' @export
metadata_table <- function(run, metadata) {
    df <- metadata[, unique(c("description",
                        "samplename",
                        run$intgroup))] %>%
        remove_rownames
    return(kable(df, caption = "Sample metadata"))
}
