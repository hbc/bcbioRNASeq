#' @rdname import
#' @description Import a project data file
#'
#' @import readr
#' @import stringr
#' @import tibble
#'
#' @param file File name
#' @param row_names Column identifier to use for row names
#' @param ... Optional parameters for \code{readr}
#'
#' @return bcbio run data
#' @export
import_file <- function(
    run,
    file,
    row_names = NULL,
    ...) {
    check_run(run)

    # Check that file exists
    filepath <- file.path(run$project_dir, file)
    if (!file.exists(filepath)) {
        stop("file could not be found")
    }

    # Detect file extension
    if (grepl("\\.[a-z]+$", file)) {
        # ext <- gsub("^.*\\.([a-z]+)$", "\\1", file)
        ext <- str_match(file, "\\.([a-z]+)$")[2]
    } else {
        stop("file does not have an extension")
    }

    # File import
    if (ext == "csv") {
        data <- read_csv(filepath, col_types = cols(), ...)
    } else if (ext == "tsv") {
        data <- read_tsv(filepath, col_types = cols(), ...)
    } else if (ext == "counts") {
        data <- read_tsv(filepath, col_types = cols(), ...) %>% as.matrix
    } else {
        stop("unsupported file extension")
    }

    # Coerce tibble to data frame. Might want to disable this in a future update
    # as tibble becomes more standardized and the rownames issue is sorted out.
    if (is_tibble(data)) {
        data <- as.data.frame(data)
    }

    # Set row names
    if (!is.null(row_names)) {
        rownames(data) <- data[[row_names]]
        data[[row_names]] <- NULL
    }

    return(data)
}
