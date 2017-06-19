#' Object integrity checks
#'
#' @rdname internal-integrity_checks
#' @keywords internal
#'
#' @param stop Stop upon check failure.
#'
#' @author Michael Steinbaugh



#' @rdname internal-integrity_checks
#'
#' @param bcb [bcbioRnaDataSet].
.check_bcb <- function(bcb) {
    if (!is.list(run)) {
        stop("Run object is not a list")
    }
    if (is.null(run$yaml)) {
        stop("Run does contain YAML info, please resave")
    }
    if (is.null(run$metadata)) {
        stop("Run does not contain metadata, please resave")
    }
    if (run$analysis == "rnaseq" & is.null(run$ensembl)) {
        message("Run does not contain Ensembl annotations")
    }

    # Metadata format and description name checks
    if (is_tibble(run$metadata) | !is.data.frame(run$metadata)) {
        stop("Run metadata not saved as a data frame")
    }
    if (!identical(rownames(run$metadata), run$metadata$description)) {
        stop("Run metadata rownames must match the description")
    }
    if (!identical(names(run$sample_dirs), run$metadata$description)) {
        stop("Run metadata descriptions don't match sample directories")
    }
    if (!any(sapply(run$metadata, is.factor))) {
        stop("Run metadata does not contain factors")
    }
}



#' @rdname internal-integrity_checks
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



#' @rdname internal-integrity_checks
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



#' @rdname internal-integrity_checks
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


#' @rdname internal-integrity_checks
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
