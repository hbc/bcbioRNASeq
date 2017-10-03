#' Read Program Versions
#'
#' @familly Project Directory File Utilities
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param projectDir Project directory.
#'
#' @return [data.frame].
.programs <- function(projectDir) {
    file <- file.path(projectDir, "programs.txt")
    if (!file.exists(file)) {
        warning(paste(basename(file), "missing"))
        return(NULL)
    }
    read_delim(file, col_names = c("program", "version"), delim = ",")
}
