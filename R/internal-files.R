#' Read file
#'
#' Supports automatic loading of `.csv`, `.tsv`, `.xlsx`, and `.counts` files.
#'
#' @rdname read_file
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param file File path.
#' @param column_to_rownames Column identifier to use for row names.
#' @param ... Additional parameters.
#'
#' @return [DataFrame].
#'
#' @seealso
#' - [readr](http://readr.tidyverse.org)
#' - [readxl](http://readxl.tidyverse.org)
.read_file <- function(file, column_to_rownames = NULL, ...) {
    if (is.null(file)) {
        return(NULL)
    }
    if (!is.character(file)) {
        stop("File path must be a string")
    }

    file_path <- normalizePath(file)
    file_name <- basename(file_path)

    if (!file.exists(file_path)) {
        stop(paste(file_name, "not found"))
    }

    message(paste("Reading", file_name))

    # Detect file extension
    if (grepl("\\.[a-z]+$", file_name)) {
        ext <- str_match(file_name, "\\.([a-z]+)$")[[2L]]
    } else {
        stop("File extension missing")
    }

    # File import, based on extension
    if (ext == "csv") {
        data <- read_csv(file_path, ...)
    } else if (ext == "tsv") {
        data <- read_tsv(file_path, ...)
    } else if (ext == "txt") {
        data <- read_delim(file_path, ...)
    } else if (ext == "xlsx") {
        data <- read_excel(file_path, ...)
    } else if (ext == "counts") {
        data <- read_tsv(file_path, ...) %>% as.matrix
    } else {
        stop("Unsupported file type")
    }

    # Coerce tibble to data frame, to allow for rownames
    if (is_tibble(data)) {
        data <- as.data.frame(data)
    }

    # Set row names, if desired
    if (!is.null(column_to_rownames)) {
        rownames(data) <- data[[column_to_rownames]]
        data[[column_to_rownames]] <- NULL
    }

    # Finally, strip all NA columns and rows, then set as S4 DataFrame
    data %>%
        remove_na %>%
        DataFrame
}



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
    meta <- .read_file(file)
    if (is.null(meta)) {
        return(NULL)
    }

    # Coerce to data.frame
    meta <- as.data.frame(meta)

    # First column must be the FASTQ file name
    names(meta)[[1L]] <- "file_name"

    meta <- meta %>%
        # Strip all NA rows and columns
        remove_na %>%
        # Make names snake_case
        snake %>%
        # Remove rows with no description
        filter(!is.na(.data[["description"]]))

    # Lane split, if desired
    if (is.numeric(lanes)) {
        meta <- meta %>%
            group_by(!!sym("file_name")) %>%
            # Expand by lane (e.g. "L001")
            expand_(~paste0("L", str_pad(1L:lanes, 3L, pad = "0"))) %>%
            left_join(meta, by = "file_name") %>%
            ungroup %>%
            mutate(file_name = paste(.data[["file_name"]],
                                     .data[["lane"]],
                                     sep = "_"),
                   description = .data[["file_name"]])
    }

    # Subset by pattern, if desired
    if (!is.null(pattern)) {
        meta <- meta[str_detect(meta[[pattern_col]], pattern), ]
    }

    # Convert to data frame, coerce to factors, and set rownames
    meta %>%
        mutate_all(factor) %>%
        mutate(file_name = as.character(.data[["file_name"]]),
               description = as.character(.data[["description"]])) %>%
        arrange(!!sym("description")) %>%
        DataFrame %>%
        set_rownames(.[["description"]])
}



#' Read data versions
#'
#' @rdname data_versions
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param project_dir Project directory.
.data_versions <- function(project_dir) {
    file.path(project_dir, "data_versions.csv") %>% .read_file
}



#' Read program versions
#'
#' @rdname program_versions
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param project_dir Project directory.
.programs <- function(project_dir) {
    file.path(project_dir, "programs.txt") %>%
        .read_file(col_names = c("program", "version"), delim = ",")
}
