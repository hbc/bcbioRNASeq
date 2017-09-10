#' Read Data Versions
#'
#' @author Michael Steinbaugh
#' @keywords internal
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



#' Read Log File
#'
#' @author Michael Steinbaugh
#' @keywords internal
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



#' Read Program Versions
#'
#' @author Michael Steinbaugh
#' @keywords internal
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
