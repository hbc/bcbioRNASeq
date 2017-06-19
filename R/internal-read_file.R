#' Read custom metadata file
#'
#' @rdname custom_metadata
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param file Metadata file. CSV and XLSX formats are supported.
#' @param pattern Apply grep pattern matching to samples
#' @param pattern_col Column in data frame used for pattern subsetting
#' @param lanes Number of lanes used to split the samples into technical
#'   replicates (`_LXXX`) suffix.
#'
#' @return [DataFrame].
.custom_metadata <- function(
    file,
    pattern = NULL,
    pattern_col = "description",
    lanes = NULL) {
    if (!file.exists(file)) {
        stop("File not found")
    }

    if (grepl("\\.xlsx$", file)) {
        metadata <- read_excel(file)
    } else {
        metadata <- read_csv(file)
    }

    # First column must be the FASTQ file name
    names(metadata)[1] <- "file_name"

    metadata <- metadata %>%
        as_tibble %>%
        # Strip all NA rows and columns
        remove_na %>%
        # Make names snake_case
        snake %>%
        # Remove rows with no description
        filter(!is.na(.data$description))

    # Lane split, if desired
    if (is.numeric(lanes)) {
        metadata <- metadata %>%
            group_by(!!sym("file_name")) %>%
            # Expand by lane (e.g. "L001")
            expand_(~paste0("L", str_pad(1:lanes, 3, pad = "0"))) %>%
            left_join(metadata, by = "file_name") %>%
            ungroup %>%
            mutate(file_name = paste(.data[["file_name"]],
                                     .data[["lane"]],
                                     sep = "_"),
                   description = .data[["file_name"]])
    }

    # Subset by pattern, if desired
    if (!is.null(pattern)) {
        metadata <- metadata[str_detect(metadata[[pattern_col]], pattern), ]
    }

    # Convert to data frame, coerce to factors, and set rownames
    metadata %>%
        mutate_all(factor) %>%
        mutate(file_name = as.character(.data$file_name),
               description = as.character(.data$description)) %>%
        arrange(!!sym("description")) %>%
        as.data.frame %>%
        set_rownames(.$description) %>%
        DataFrame
}



#' Read file from bcbio-nextgen final upload directory
#'
#' @rdname file
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param parent_dir Parent directory.
#' @param file_name File name.
#' @param column_to_rownames Column identifier to use for row names.
#' @param ... Passthrough parameters for [read_csv()] or [read_tsv()].
#'
#' @return Miscellaneous data frame.
#' @export
.file <- function(
    parent_dir,
    file_name,
    column_to_rownames = NULL,
    ...) {
    # Check that file exists
    file_path <- file.path(project_dir, file_name)
    if (!file.exists(file_path)) {
        return(NULL)
    }

    # Detect file extension
    if (grepl("\\.[a-z]+$", file_name)) {
        # ext <- gsub("^.*\\.([a-z]+)$", "\\1", file_name)
        ext <- str_match(file_name, "\\.([a-z]+)$")[2]
    } else {
        stop("File does not have an extension")
    }

    # File import
    if (ext == "csv") {
        data <- read_csv(file_path, ...)
    } else if (ext == "tsv") {
        data <- read_tsv(file_path, ...)
    } else if (ext == "counts") {
        data <- read_tsv(file_path, ...) %>% as.matrix
    } else {
        stop("Unsupported file extension")
    }

    # Coerce tibble to data frame
    if (is_tibble(data)) {
        data <- as.data.frame(data)
    }

    # Set row names, if desired
    if (!is.null(column_to_rownames)) {
        rownames(data) <- data[[column_to_rownames]]
        data[[column_to_rownames]] <- NULL
    }

    data %>% remove_na
}
