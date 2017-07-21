#' Read Data Versions
#'
#' @rdname data_versions
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param project_dir Project directory.
.data_versions <- function(project_dir) {
    file <- file.path(project_dir, "data_versions.csv")
    if (!file.exists(file)) {
        warning("Data versions file missing")
        return(NULL)
    }
    read_csv(file)
}



#' Read Log File
#'
#' @rdname log_file
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param file Log file.
.log_file <- function(file) {
    if (!file.exists(file)) return(NULL)
    read_lines(file)
}



#' Read Program Versions
#'
#' @rdname program_versions
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param project_dir Project directory.
.programs <- function(project_dir) {
    file <- file.path(project_dir, "programs.txt")
    if (!file.exists(file)) {
        warning("Program version file missing")
        return(NULL)
    }
    read_delim(file, col_names = c("program", "version"), delim = ",")
}
