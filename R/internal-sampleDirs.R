#' Detect Sample Directories
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom fs dir_ls path_real
#'
#' @param uploadDir Upload directory.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will abort if no sample directories match.
#' @noRd
.sampleDirs <- function(uploadDir) {
    assert_all_are_dirs(uploadDir)
    subdirs <- dir_ls(
        path = uploadDir,
        type = "directory",
        recursive = FALSE
    )
    subdirPattern <- paste0(perSampleDirs, collapse = "|") %>%
        paste0("^", ., "$")
    sampleDirs <- dir_ls(
        path = subdirs,
        recursive = FALSE,
        regexp = subdirPattern
    )
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
    sampleDirs <- path_real(sampleDirs)
    names(sampleDirs) <- names

    inform(paste(length(sampleDirs), "samples detected"))

    sampleDirs
}
