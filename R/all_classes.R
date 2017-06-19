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
bcbioRnaDataSet <- setClass(
    "bcbioRnaDataSet", contains = "SummarizedExperiment",
    representation = representation(callers="SimpleList"))
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
    bcb <- new("bcbioRnaDataSet", se)
    bcbio(bcb, "tximport") <- run$txi
    metadata(bcb)[["txi"]] <- NULL
    # [! Fix] copy real data to new pointer so we can remove it from here
    # run$txi <- NULL
    bcbio(bcb, "featureCounts") <- read_bcbio_file(metadata(bcb), "combined.counts")
    if (run$analysis == "srnaseq"){
        bcb <- read_smallrna_counts(bcb)
    }
    bcb
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
        intgroup = intgroup,
        metadata_file = metadata_file,
        organism = organism)
    genes <- as.data.frame(run$txi$counts) %>%
        rownames_to_column("ensembl_gene_id") %>%
        left_join(run$ensembl)
    row.names(genes) <- row.names(run$txi$counts)
    se <- SummarizedExperiment(
        assays = SimpleList(counts = run$txi$counts),
        rowData = DataFrame(genes),
        colData = DataFrame(run$metadata),
        metadata = run, ...)
    .bcbioRnaDataSet(se, run)
}
