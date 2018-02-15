#' Detect Sample Directories
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param uploadDir Upload directory.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will [abort()] if no complete sample directories match.
#' @noRd
.sampleDirs <- function(uploadDir) {
    assert_all_are_dirs(uploadDir)
    subdirs <- list.dirs(uploadDir, full.names = TRUE, recursive = FALSE)
    subdirPattern <- paste0(perSampleDirs, collapse = "|") %>%
        paste0("^", ., "$")
    sampleDirs <- list.files(
        subdirs,
        pattern = subdirPattern,
        full.names = TRUE,
        recursive = FALSE)
    assert_is_non_empty(sampleDirs)
    sampleDirs <- sampleDirs %>%
        dirname() %>%
        sort() %>%
        unique()

    # Ensure removal of nested `projectDir`
    if (any(grepl(projectDirPattern, basename(sampleDirs)))) {
        sampleDirs <- sampleDirs %>%
            .[!grepl(projectDirPattern, basename(.))]
        assert_is_non_empty(sampleDirs)
    }

    # Generate names from file paths and make valid
    names <- basename(sampleDirs) %>%
        make.names(unique = TRUE) %>%
        gsub("\\.", "_", .)
    sampleDirs <- normalizePath(sampleDirs)
    names(sampleDirs) <- names

    inform(paste(length(sampleDirs), "samples detected"))

    sampleDirs
}
