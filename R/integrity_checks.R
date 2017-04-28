#' Perform object integrity checks
#'
#' @rdname integrity_checks
#' @author Michael Steinbaugh



#' @rdname integrity_checks
#' @description Check \code{bcbio-nextgen} run
#' @param run \code{bcbio-nextgen} run
#' @export
check_run <- function(run) {
    if (!is.list(run)) {
        stop("Run object is not a list")
    }
    if (!length(dir(run$final_dir))) {
        stop("Could not access run_dir")
    }
    if (!length(run$sample_dirs)) {
        stop("No sample directories in run")
    }
    if (!length(dir(run$project_dir))) {
        stop("Could not access project_dir")
    }
    if (is.null(run$metadata)) {
        stop("Run does not contain metadata, please reload")
    }
    if (is.null(run$yaml)) {
        stop("Run does contain YAML info, please reload")
    }
    # Ensure metadata descriptions match the sample folders
    if (!identical(run$metadata$description, names(run$sample_dirs))) {
        stop("Metadata descriptions don't match the sample names")
    }
}



#' @rdname integrity_checks
#' @description Check \code{DESeqDataSet}
#' @param dds \code{DESeqDataSet}
#' @param stop Stop if class doesn't match
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
#' @description Check \code{DESeqTransform}
#' @param dt \code{DESeqTransform}
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
#' @description Check \code{DESeqResults}
#' @param res \code{DESeqResults}
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
