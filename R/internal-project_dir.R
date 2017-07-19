#' Read Data Versions
#'
#' @rdname data_versions
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param project_dir Project directory.
.data_versions <- function(project_dir) {
    file.path(project_dir, "data_versions.csv") %>% read_csv
}



#' Read Program Versions
#'
#' @rdname program_versions
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param project_dir Project directory.
.programs <- function(project_dir) {
    file.path(project_dir, "programs.txt") %>%
        read_delim(col_names = c("program", "version"),
                   delim = ",")
}
