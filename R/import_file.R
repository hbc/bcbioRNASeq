#' Import a bcbio project data file
#'
#' @author Michael Steinbaugh
#'
#' @import readr
#' @import stringr
#' @import tibble
#'
#' @param bcbio bcbio run object
#' @param file File name
#' @param rownames Column identifier to use for rownames
#' @param ... Optional pass through to \code{readr} package
#'
#' @return bcbio data
#' @export
#'
#' @examples
#' \dontrun{
#' import_file(bcbio, file = "tx2gene.csv", rownames = "id")
#' }
import_file <- function(bcbio,
                        file,
                        rownames = NULL) {
    filepath <- file.path(bcbio$project_dir, file)
    if (!file.exists(filepath)) {
        stop("File could not be found.")
    }

    # Detect file extension
    if (grepl("\\.[a-z]+$", file)) {
        # base method
        # ext <- gsub("^.*\\.([a-z]+)$", "\\1", file)
        ext <- stringr::str_match(file, "\\.([a-z]+)$")[2]
    } else {
        stop("File does not have an extension.")
    }

    # File import
    if (ext == "csv") {
        data <- readr::read_csv(filepath, ...)
    } else if (ext == "tsv") {
        data <- readr::read_tsv(filepath, ...)
    } else {
        stop("Unsupported file extension.")
    }

    # Coerce to data frame. Might want to disable this in a future update as
    # tibble becomes more standardized and the rownames issue is sorted out.
    if (tibble::is_tibble(data)) {
        data <- as.data.frame(data)
    }

    # Set rownames
    if (!is.null(rownames)) {
        rownames(data) <- data[[rownames]]
        data[[rownames]] <- NULL
    }

    return(data)
}
