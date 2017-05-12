#' Perform object integrity checks.
#'
#' @rdname integrity_checks
#' @author Michael Steinbaugh



#' @rdname integrity_checks
#' @keywords internal
#' @description Check bcbio-nextgen run. Always stop on failure.
#' @param run bcbio-nextgen run
#' @export
check_run <- function(run) {
    if (!is.list(run)) {
        stop("Run object is not a list")
    }
    if (is.null(run$yaml)) {
        stop("Run does contain YAML info, please resave")
    }
    if (!length(dir(run$upload_dir))) {
        stop("Could not access upload dir")
    }
    if (!length(dir(run$project_dir))) {
        stop("Could not access project dir")
    }
    if (!length(run$sample_dirs)) {
        stop("No sample directories in run")
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
#' @keywords internal
#' @description Check \linkS4class{DESeqDataSet}.
#' @param dds \linkS4class{DESeqDataSet}.
#' @param stop Stop if class doesn't match.
#' @export
check_dds <- function(dds, stop = TRUE) {
    if (class(dds)[1] == "DESeqDataSet") {
        return(TRUE)
    } else {
        if (isTRUE(stop)) {
            stop("DESeqDataSet required")
        } else {
            return(FALSE)
        }
    }
}



#' @rdname integrity_checks
#' @description Check \linkS4class{DESeqTransform}.
#' @param dt \linkS4class{DESeqTransform}.
#' @export
check_dt <- function(dt, stop = TRUE) {
    if (class(dt)[1] == "DESeqTransform") {
        return(TRUE)
    } else {
        if (isTRUE(stop)) {
            stop("DESeqDataSet required")
        } else {
            return(FALSE)
        }
    }
}



#' @rdname integrity_checks
#' @description Check \linkS4class{DESeqResults}.
#' @param res \linkS4class{DESeqResults}.
#' @export
check_res <- function(res, stop = TRUE) {
    if (class(res)[1] == "DESeqResults") {
        return(TRUE)
    } else {
        if (isTRUE(stop)) {
            stop("DESeqResults required")
        } else {
            return(FALSE)
        }
    }
}
