#' Read metadata
#'
#' @param csv Comma separated values file
#' @param pattern Apply grep pattern matching to samples
#' @param pattern_col Column in data frame used for pattern subsetting
#'
#' @return Metadata data frame
#' @export
read_metadata <- function(
    csv,
    pattern = NULL,
    pattern_col = "description") {
    metadata <- read_csv(csv, col_types = cols()) %>%
        set_names_snake %>%
        arrange(!!sym("description"))
    # Subset by pattern, if desired
    if (!is.null(pattern)) {
        metadata <- metadata[str_detect(metadata[[pattern_col]], pattern), ]
    }
    return(metadata)
}
