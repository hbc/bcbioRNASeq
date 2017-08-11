#' Read Custom Metadata File
#'
#' @rdname internal-readSampleMetadataFile
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param file Metadata file. CSV and XLSX formats are supported.
#' @param pattern *Optional*. Grep pattern to match against sample names.
#' @param patternCol *Optional*. Column in data frame used for pattern
#'   subsetting.
#' @param lanes *Optional*. Number of lanes used to split the samples into
#'   technical replicates (`_LXXX`) suffix.
#'
#' @return [tibble] grouped by `sampleName`.
.readSampleMetadataFile <- function(
    file,
    pattern = NULL,
    patternCol = "sampleName",
    lanes = 1L) {
    if (is.null(file)) return(NULL)
    meta <- readFileByExtension(file) %>% as("tibble")
    # First column must match the sample directory names output by bcbio.
    # This is set as `description` during bcbio config. Any other column
    # is optional metadata.
    names(meta)[[1L]] <- "sampleName"
    meta <- meta %>%
        # Strip all NA rows and columns
        removeNA %>%
        # Make colnames camelCase
        camel %>%
        # Remove rows with no sample name. Sometimes Excel files will add
        # empty rows, so this helps correct that problem as well.
        .[!is.na(.[["sampleName"]]), ]

    # Lane split, if desired
    if (lanes > 1L) {
        meta <- meta %>%
            group_by(!!sym("sampleName")) %>%
            # Expand by lane (e.g. "L001")
            expand_(dots = ~paste0("L", str_pad(1L:lanes, 3L, pad = "0"))) %>%
            # `expand_cols` param doesn't seem to work in tidyr 0.6.3, so
            # set manually here instead
            set_colnames(c("sampleName", "lane")) %>%
            left_join(meta, by = "sampleName") %>%
            ungroup %>%
            mutate(sampleName = paste(.data[["sampleName"]],
                                      .data[["lane"]],
                                      sep = "_"),
                   sampleID = .data[["sampleName"]])
    }

    # Subset by pattern, if desired
    if (!is.null(pattern)) {
        meta <- meta %>%
            .[str_detect(.[[patternCol]], pattern), ]
    }

    meta %>%
        .metaPriorityCols %>%
        .metaFactors
}
