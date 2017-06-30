#' Read file
#'
#' Supports automatic loading of standard `.csv`, `.tsv`, and `.xlsx`
#' files. Also supports bcbio-nextgen pipeline-specific `.counts`, `.colnames`,
#' and `.rownames` files.
#'
#' @rdname read_file
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param file File path.
#' @param ... Additional parameters.
#'
#' @return [tibble].
#'
#' @seealso [readr](http://readr.tidyverse.org),
#'   [readxl](http://readxl.tidyverse.org).
.read_file <- function(file, ...) {
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
        data <- read_csv(file_path, progress = FALSE, ...)
    } else if (ext == "tsv") {
        data <- read_tsv(file_path, progress = FALSE, ...)
    } else if (ext == "txt") {
        data <- read_delim(file_path, progress = FALSE, ...)
    } else if (ext == "xlsx") {
        data <- read_excel(file_path, ...)
    } else if (ext %in% c("colnames", "rownames")) {
        data <- read_lines(file_path, ...)
    } else if (ext == "counts") {
        data <- read_tsv(file_path, progress = FALSE, ...)
    } else {
        stop("Unsupported file type")
    }

    # Coerce data frame to tibble, if necessary
    if (is.data.frame(data) & !is_tibble(data)) {
        data <- as_tibble(data)
    }

    # Return
    if (is_tibble(data)) {
        data %>%
            remove_na %>%
            snake(rownames = FALSE)
    } else {
        data
    }
}



#' Read data versions
#'
#' @rdname data_versions
#' @keywords internal
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
#' @author Michael Steinbaugh
#'
#' @param project_dir Project directory.
.programs <- function(project_dir) {
    file.path(project_dir, "programs.txt") %>%
        .read_file(col_names = c("program", "version"), delim = ",")
}
