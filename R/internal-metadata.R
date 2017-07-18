.meta_priority_cols <- function(meta) {
    meta %>%
        as("tibble") %>%
        # Sanitize `sample_id` into snake_case
        mutate(sample_id = snake(.data[["sample_id"]])) %>%
        tidy_select(unique(c(meta_priority_cols, sort(colnames(.))))) %>%
        arrange(!!!syms(meta_priority_cols))
}

.meta_factors <- function(meta) {
    meta %>%
        mutate_if(!colnames(.) %in% meta_priority_cols, factor)
}



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
#' @return [tibble] grouped by `sample_name`.
.sample_metadata_file <- function(
    file,
    pattern = NULL,
    pattern_col = "sample_name",
    lanes = 1L) {
    if (is.null(file)) return(NULL)
    meta <- read_file_by_extension(file) %>% as("tibble")
    # First column must be the `sample_name`, which points to the FASTQ.
    # bcbio labels this `samplename` by default. Rename to `sample_name`
    # here to ensure consistent snake_case naming syntax.
    names(meta)[[1L]] <- "sample_id"
    meta <- meta %>%
        # Strip all NA rows and columns
        remove_na %>%
        # Make colnames snake_case
        snake %>%
        # Remove rows with no description. Sometimes Excel files will add
        # empty rows, so this helps correct that problem as well.
        filter(!is.na(.data[["sample_name"]]))

    # Lane split, if desired
    if (lanes > 1) {
        meta <- meta %>%
            group_by(!!sym("sample_name")) %>%
            # Expand by lane (e.g. "L001")
            expand_(dots = ~str_c("L", str_pad(1L:lanes, 3L, pad = "0"))) %>%
            # `expand_cols` param doesn't seem to work in tidyr 0.6.3, so
            # set manually here instead
            set_colnames(c("sample_name", "lane")) %>%
            left_join(meta, by = "sample_name") %>%
            ungroup %>%
            mutate(sample_name = paste(.data[["sample_name"]],
                                     .data[["lane"]],
                                     sep = "_"),
                   sample_id = .data[["sample_name"]])
    }

    # Subset by pattern, if desired
    if (!is.null(pattern)) {
        meta <- meta %>%
            filter(str_detect(.data[[pattern_col]], pattern))
    }

    meta %>% .meta_priority_cols %>% .meta_factors
}
