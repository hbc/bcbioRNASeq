#' Read data versions
#' @rdname internal-data_versions
#' @param project_dir Project directory.
.data_versions <- function(project_dir) {
    file <- file.path(project_dir, "data_versions.csv")
    if (!file.exists(file)) {
        return(NULL)
    }
    message("Reading data versions")
    read_csv(file,
             # c = character
             # T = date time
             col_types = "ccT") %>% DataFrame
}



#' Read program versions
#' @rdname internal-program_versions
#' @param project_dir Project directory.
.program_versions <- function(project_dir) {
    file <- file.path(project_dir, "programs.txt")
    if (!file.exists(file)) {
        return(NULL)
    }
    message("Reading program versions")
    read_delim(file,
               delim = ",",
               col_names = c("program", "version"),
               # c = character
               col_types = "cc") %>% DataFrame
}
