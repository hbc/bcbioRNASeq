#' Read Data Versions
#'
#' @familly Project Directory File Utilities
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param projectDir Project directory.
#'
#' @return [data.frame].
.dataVersions <- function(projectDir) {
    file <- file.path(projectDir, "data_versions.csv")
    if (!file.exists(file)) {
        warning(paste(basename(file), "missing"))
        return(NULL)
    }
    read_csv(file)
}
