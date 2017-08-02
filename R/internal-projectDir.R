#' Read Data Versions
#'
#' @rdname internal-dataVersions
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param projectDir Project directory.
.dataVersions <- function(projectDir) {
    file <- file.path(projectDir, "data_versions.csv")
    if (!file.exists(file)) {
        warning(paste(basename(file), "missing"))
        return(NULL)
    }
    read_csv(file)
}



#' Read Log File
#'
#' @rdname internal-logFile
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param file Log file.
.logFile <- function(file) {
    if (!file.exists(file)) {
        warning(paste(basename(file), "missing"))
        return(NULL)
    }
    read_lines(file)
}



#' Read Program Versions
#'
#' @rdname internal-programVersions
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param projectDir Project directory.
.programs <- function(projectDir) {
    file <- file.path(projectDir, "programs.txt")
    if (!file.exists(file)) {
        warning(paste(basename(file), "missing"))
        return(NULL)
    }
    read_delim(file, col_names = c("program", "version"), delim = ",")
}
