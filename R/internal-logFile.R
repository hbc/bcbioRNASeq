#' Read Log File
#'
#' @family Project Directory File Utilities
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param file Log file.
#'
#' @return Character vector.
.logFile <- function(file) {
    if (!file.exists(file)) {
        warning(paste(basename(file), "missing"))
        return(NULL)
    }
    read_lines(file)
}
