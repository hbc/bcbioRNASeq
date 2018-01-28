#' Detect Sample Directories
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom stats setNames
#'
#' @param uploadDir Upload directory.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will [abort()] if no complete sample directories match.
#' @noRd
.sampleDirs <- function(uploadDir) {
    subdirs <- list.dirs(uploadDir, full.names = TRUE, recursive = FALSE)
    subdirPattern <- paste0(perSampleDirs, collapse = "|") %>%
        paste0("^", ., "$")
    sampleDirs <- list.files(
        subdirs,
        pattern = subdirPattern,
        full.names = TRUE,
        recursive = FALSE) %>%
        dirname() %>%
        sort() %>%
        unique()

    # Ensure removal of nested `projectDir`
    if (any(grepl(x = basename(sampleDirs), pattern = projectDirPattern))) {
        sampleDirs <- sampleDirs %>%
            .[!grepl(x = basename(.), pattern = projectDirPattern)]
    }

    # Return
    if (length(sampleDirs) == 0L) {
        abort("No sample directories detected")
    } else {
        # Generate names from file paths and make valid
        names <- basename(sampleDirs) %>%
            make.names(unique = TRUE) %>%
            gsub(x = ., pattern = "\\.", replacement = "_")
        sampleDirs <- normalizePath(sampleDirs) %>%
            setNames(names)
        inform(paste(length(sampleDirs), "samples detected"))
    }
    sampleDirs
}
