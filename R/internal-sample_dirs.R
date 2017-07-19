#' Detect Sample Directories
#'
#' @rdname sample_dirs
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param upload_dir Upload directory.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will [stop()] if no complete sample directories match.
.sample_dirs <- function(upload_dir) {
    subdirs <- list.dirs(upload_dir, full.names = TRUE, recursive = FALSE)
    subdir_pattern <- str_c(per_sample_dirs, collapse = "|") %>%
        str_c("^", ., "$")
    sample_dirs <- list.files(subdirs,
               pattern = subdir_pattern,
               full.names = TRUE,
               recursive = FALSE) %>%
        dirname %>%
        sort %>%
        unique

    # Ensure removal of nested `project_dir`
    if (any(str_detect(basename(sample_dirs), project_dir_pattern))) {
        sample_dirs <- sample_dirs %>%
            .[!str_detect(basename(.), project_dir_pattern)]
    }

    # Return
    if (length(sample_dirs) == 0L) {
        stop("No sample directories detected")
    } else {
        # Generate names from file paths and sanitize into snake_case
        names <- basename(sample_dirs) %>% snake
        sample_dirs <- normalizePath(sample_dirs) %>%
            set_names(names)
        message(paste(length(sample_dirs), "samples detected"))
    }
    sample_dirs
}
