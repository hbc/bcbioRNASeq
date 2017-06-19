#' Object integrity checks
#'
#' @rdname integrity
#' @keywords internal
#'
#' @param stop Stop upon check failure.
#'
#' @author Michael Steinbaugh



#' @rdname integrity
#' @param dds [DESeqDataSet].
.check_dds <- function(dds, stop = TRUE) {
    if (class(dds)[1] == "DESeqDataSet") {
        TRUE
    } else {
        if (isTRUE(stop)) {
            stop("DESeqDataSet required")
        } else {
            FALSE
        }
    }
}



#' @rdname integrity
#' @param dt [DESeqTransform].
.check_dt <- function(dt, stop = TRUE) {
    if (class(dt)[1] == "DESeqTransform") {
        TRUE
    } else {
        if (isTRUE(stop)) {
            stop("DESeqDataSet required")
        } else {
            FALSE
        }
    }
}



#' @rdname integrity
#' @param res [DESeqResults].
.check_res <- function(res, stop = TRUE) {
    if (class(res)[1] == "DESeqResults") {
        TRUE
    } else {
        if (isTRUE(stop)) {
            stop("DESeqResults required")
        } else {
            FALSE
        }
    }
}



#' @rdname integrity
#' @param sample_dirs Named character vector of sample directory paths.
.check_sample_dirs <- function(sample_dirs) {
    if (!is.character(sample_dirs)) {
        stop("sample_dirs must be a character vector")
    }
    # Check for named character vector
    if (is.null(names(sample_dirs))) {
        stop("sample_dirs must have assigned sample names")
    }
    if (!identical(basename(sample_dirs), sample_names)) {
        stop("Sample name assignment mismatch")
    }
}
