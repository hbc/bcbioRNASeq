#' Perform object integrity checks
#'
#' @rdname integrity_checks
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param stop Stop upon check failure.



#' @rdname integrity_checks
#' @param run bcbio-nextgen run.
#' @export
check_run <- function(run) {
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
        stop("RNA-seq run does not contain Ensembl annotations, please resave")
    }

    # Metadata format and description name checks
    if (is_tibble(run$metadata) | !is.data.frame(run$metadata)) {
        stop("Metadata must be saved as a data frame")
    }
    if (!identical(rownames(run$metadata), run$metadata$description)) {
        stop("Metadata rownames must match the description")
    }
    if (!identical(run$metadata$description, names(run$sample_dirs))) {
        stop("Metadata descriptions don't match the run sample directories")
    }
}



#' @rdname integrity_checks
#' @param dds [DESeqDataSet].
#' @export
check_dds <- function(dds, stop = TRUE) {
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



#' @rdname integrity_checks
#' @param dt [DESeqTransform].
#' @export
check_dt <- function(dt, stop = TRUE) {
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



#' @rdname integrity_checks
#' @param res [DESeqResults].
#' @export
check_res <- function(res, stop = TRUE) {
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
