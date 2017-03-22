#' Import a bcbio project data file
#'
#' @author Michael Steinbaugh
#'
#' @import readr
#'
#' @param bcbio bcbio run object
#' @param filename Filename
#' @param input Input format
#' @param output Desired output format
#' @param rownames Column identifier to use for rownames
#'
#' @return bcbio data
#' @export
#'
#' @examples
#' \dontrun{
#' import_file(bcbio, "combined.counts",
#'             output = "matrix", rownames = "id")
#' }
import_file <- function(bcbio,
                        filename,
                        input = "tsv",
                        output = "data.frame",
                        rownames = NULL) {
    filepath <- file.path(bcbio$project_dir, filename)
    if (!file.exists(filepath)) {
        stop("File could not be found.")
    }

    # File import
    if (input == "csv") {
        data <- readr::read_csv(filepath)
    } else if (input == "tsv") {
        data <- readr::read_tsv(filepath)
    }

    # Coerce to data frame
    data <- as.data.frame(data)

    # Set rownames
    if (!is.null(rownames)) {
        rownames(data) <- data[[rownames]]
        data[[rownames]] <- NULL
    }

    # Convert to desired output
    if (output == "matrix") {
        data <- as.matrix(data)
    }

    return(data)
}
