#' Class that contains bcbio run information
#'
#' `bcbioRnaDataSet` is a subclass of [SummarizedExperiment] designed to store
#' an RNA-seq analysis. This class contains experiment metadata, raw counts,
#' normalilzed counts, and summary statistics for each sample analyzed.
#'
#' Methods for this objects ...
#'
#' `metadata` contains ...
#'
#' @rdname bcbioRnaDataSet
#' @keywords internal
#' @aliases bcbioRnaDataSet-class
#' @export
bcbioRnaDataSet <- setClass("bcbioRnaDataSet",
                            contains = "SummarizedExperiment")
setValidity("bcbioRnaDataSet", function(object) { TRUE })



# Constructor
.bcbioRnaDataSet <- function(se, run) {
    if (!is(se, "SummarizedExperiment")) {
        if (is(se, "SummarizedExperiment0")) {
            se <- as(se, "SummarizedExperiment")
        } else if (is(se, "SummarizedExperiment")) {
            # Only to help transition from SummarizedExperiment to new
            # RangedSummarizedExperiment objects, remove once transition is
            # complete
            se <- as(se, "SummarizedExperiment")
        } else {
            stop("'se' must be a SummarizedExperiment object")
        }
    }
    ids <- new("bcbioRnaDataSet", se)
    metadata(ids) <- run
    ids
}



#' @rdname load_run
#' @keywords internal
#' @param ... Arguments provided to [SummarizedExperiment]. Remove this once
#' we migrate the S4 loading variant over as main method.
#' @export
load_run_S4 <- function(
    upload_dir = "final",
    analysis = "rnaseq",
    intgroup = "description",
    metadata_file = NULL,
    organism = NULL,
    ...) {
    run <- load_run(
        upload_dir,
        analysis,
        intgroup,
        metadata_file,
        organism)
    se <- SummarizedExperiment(
        assays = SimpleList(counts = run$txi),
        colData = DataFrame(run$metadata), ...)
    run$txi = NULL
    .bcbioRnaDataSet(se, run)
}
