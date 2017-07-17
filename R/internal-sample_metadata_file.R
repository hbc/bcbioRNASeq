#' Read custom metadata file
#'
#' @rdname sample_metadata_file
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param file Metadata file. CSV and XLSX formats are supported.
#' @param pattern *Optional*. Grep pattern to match against sample names.
#' @param pattern_col *Optional*. Column in data frame used for pattern
#'   subsetting.
#' @param lanes *Optional*. Number of lanes used to split the samples into
#'   technical replicates (`_LXXX`) suffix.
#'
#' @return [tibble] grouped by description.
.sample_metadata_file <- function(
    file,
    pattern = NULL,
    pattern_col = "description",
    lanes = NULL) {

    if (is.null(file))
        return(NULL)

    meta <- read_file_by_extension(file) %>% as("tibble")
    # First column must be the `sample_name`, which points to the FASTQ.
    # bcbio labels this `samplename` by default. Rename to `sample_name`
    # here to ensure consistent snake_case naming syntax.
    names(meta)[[1L]] <- "sample_name"
    meta <- meta %>%
        # Strip all NA rows and columns
        remove_na %>%
        # Make names snake_case
        snake %>%
        # Remove rows with no description. Sometimes Excel files will add
        # empty rows, so this helps correct that problem as well.
        filter(!is.na(.data[["description"]]))

    # Lane split, if desired
    if (is.numeric(lanes)) {
        meta <- meta %>%
            group_by(!!sym("sample_name")) %>%
            # Expand by lane (e.g. "L001")
            expand_(~paste0("L", str_pad(1L:lanes, 3L, pad = "0"))) %>%
            left_join(meta, by = "sample_name") %>%
            ungroup %>%
            mutate(sample_name = paste(.data[["sample_name"]],
                                     .data[["lane"]],
                                     sep = "_"),
                   description = .data[["sample_name"]])
    }

    # Subset by pattern, if desired
    if (!is.null(pattern)) {
        meta <- meta %>%
            filter(str_detect(.data[[pattern_col]], pattern))
    }

    # Convert to data frame, coerce to factors, and set rownames
    meta %>%
        mutate_all(factor) %>%
        mutate(sample_name = as.character(.data[["sample_name"]]),
               description = as.character(.data[["description"]])) %>%
        arrange(!!!syms(meta_priority_cols)) %>%
        group_by(!!sym("description"))
}
