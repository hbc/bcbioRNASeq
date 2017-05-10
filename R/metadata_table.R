#' Metadata table
#'
#' Returns a subset of metadata columns of interest used for knit reports. These
#' "interesting group" columns are defined as \code{intgroup} in the
#' bcbio-nextgen run object.
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run
#'
#' @return Data frame containing only the columns of interest
#' @export
metadata_table <- function(run) {
    df <- run$metadata %>% remove_rownames
    return(kable(df, caption = "Sample metadata"))
}
