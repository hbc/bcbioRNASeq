#' Read metadata
#'
#' @param csv Comma separated values file
#' @param pattern Apply grep pattern matching to samples
#' @param pattern_col Column in data frame used for pattern subsetting
#' @param save Save csv and rda files
#'
#' @return Metadata data frame
#' @export
read_metadata <- function(
    csv,
    pattern = NULL,
    pattern_col = "description",
    save = FALSE) {
    metadata <- csv %>%
        read_csv(col_types = cols()) %>%
        set_names_snake %>%
        arrange(!!sym("description"))

    # Subset by pattern, if desired
    if (!is.null(pattern)) {
        metadata <- metadata[str_detect(metadata[[pattern_col]], pattern), ]
    }

    # Save files to disk, if desired
    if (isTRUE(save)) {
        save(metadata, file = file.path("data", "metadata.rda"))
        write_csv(metadata, file.path("meta", "metadata.csv"))
    }

    return(metadata)
}
