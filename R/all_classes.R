#' Class that contains bcbio run information
#'
#' The \code{\link{bcbioRnaDataSet}} is a subclass of
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' used to store the raw data, intermediate calculations and results of an
#' miRNA/isomiR analysis. This class stores all raw isomiRs
#' data for each sample, processed information,
#' summary for each isomiR type,
#' raw counts, normalized counts, and table with
#' experimental information for each sample.
#'
#' \code{\link{bcbioRnaDataSetFromFolder}} creates this object using bcbio output files.
#'
#' Methods for this objects ...
#'
#' \code{metadata} contains ...
#'
#' @aliases bcbioRnaDataSet-class
#' @rdname bcbioRnaDataSet
#' @export
bcbioRnaDataSet <- setClass("bcbioRnaDataSet",
                          contains = "SummarizedExperiment")

setValidity( "bcbioRnaDataSet", function( object ) {
    TRUE
} )

# Constructor
.bcbioRnaDataSet <- function(se, run){
    if (!is(se, "SummarizedExperiment")) {
        if (is(se, "SummarizedExperiment0")) {
            se <- as(se, "SummarizedExperiment")
        } else if (is(se, "SummarizedExperiment")) {
            # only to help transition from SummarizedExperiment to new
            # RangedSummarizedExperiment objects,
            # remove once transition is complete
            se <- as(se, "SummarizedExperiment")
        } else {
            stop("'se' must be a SummarizedExperiment object")
        }
    }
    ids <- new("bcbioRnaDataSet", se)
    metadata(ids) <- run
    ids
}


#' \code{bcbioRnaDataSetFromFolder} loads bcbio-nextgen output
#'
#' This function parses
#' output of bcbio-nextgen tool to allow (small) RNAseq analysis of samples
#' in different groups such as
#' characterization, differential expression and clustering. It creates an
#' \code{\link[bcbioRnaseq]{bcbioRnaDataSetFromFolder}} object.
#'
#' @author Michael Steinbaugh
#' @author Lorena Pantano
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running \code{bcbio_nextgen -w template}.
#' @param analysis Analysis type (e.g. \code{rnaseq} or \code{srnaseq})
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#' @param metadata_file Custom metadata file to import. Otherwise defaults to
#'   sample metadata saved in the YAML file.
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes and normally does not need to be set.
#' @param read_counts Automatically read in the count data using
#'   \code{read_bcbio_counts()}
#'
#' @return bcbio-nextgen run object
#' @export
bcbioRnaDataSetFromFolder <- function(
    upload_dir = "final",
    analysis = "rnaseq",
    intgroup = "description",
    metadata_file = NULL,
    organism = NULL,
    read_counts = TRUE, ...) {

    run <- load_run(upload_dir, analysis, intgroup, metadata_file, organism, read_counts)
    se <- SummarizedExperiment(assays = SimpleList(counts=run$txi),
                               colData = DataFrame(run$metadata), ...)
    run$txi = NULL
    ids <- .bcbioRnaDataSet(se, run)
    return(ids)
}
