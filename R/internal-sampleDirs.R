#' Detect Sample Directories
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param uploadDir Upload directory.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will [stop()] if no complete sample directories match.
#' @noRd
.sampleDirs <- function(uploadDir) {
    subdirs <- list.dirs(uploadDir, full.names = TRUE, recursive = FALSE)
    subdirPattern <- paste0(perSampleDirs, collapse = "|") %>%
        paste0("^", ., "$")
    sampleDirs <- list.files(subdirs,
                             pattern = subdirPattern,
                             full.names = TRUE,
                             recursive = FALSE) %>%
        dirname() %>%
        sort() %>%
        unique()

    # Ensure removal of nested `projectDir`
    if (any(str_detect(basename(sampleDirs), projectDirPattern))) {
        sampleDirs <- sampleDirs %>%
            .[!str_detect(basename(.), projectDirPattern)]
    }

    # Return
    if (length(sampleDirs) == 0) {
        stop("No sample directories detected")
    } else {
        # Generate names from file paths and make valid
        names <- basename(sampleDirs) %>%
            make.names(unique = TRUE) %>%
            gsub(x = ., pattern = "\\.", replacement = "_") %>%
        sampleDirs <- normalizePath(sampleDirs) %>%
            setNames(names)
        message(paste(length(sampleDirs), "samples detected"))
    }
    sampleDirs
}
