#' Perform integrity checks on a bcbio run object
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param bcbio bcbio run object
#'
#' @return Silent, otherwise stop execution if integrity check fails
#' @export
check_bcbio_object <- function(bcbio) {
    if (!is.list(bcbio)) {
        stop("bcbio run list object not found")
    }
    if (!length(dir(bcbio$parent_dir))) {
        stop("could not access parent_dir")
    }
    if (!length(dir(bcbio$run_dir))) {
        stop("could not access run_dir")
    }
    if (!length(dir(bcbio$final_dir))) {
        stop("could not access final_dir")
    }
    if (!length(dir(bcbio$project_dir))) {
        stop("could not access project_dir")
    }
}
