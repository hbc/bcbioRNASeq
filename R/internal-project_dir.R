#' Read data versions
#'
#' @rdname data_versions
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param project_dir Project directory.
.data_versions <- function(project_dir) {
    file.path(project_dir, "data_versions.csv") %>% read_csv
}



#' Read program versions
#'
#' @rdname program_versions
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param project_dir Project directory.
.programs <- function(project_dir) {
    file.path(project_dir, "programs.txt") %>%
        read_delim(col_names = c("program", "version"),
                   delim = ",")
}
