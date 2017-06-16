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
    representation = representation(callers = "SimpleList"))
setValidity("bcbioRnaDataSet", function(object) TRUE)



# Constructor
.bcbioRnaDataSet <- function(se, run) {
    if (!is(se, "SummarizedExperiment")) {
        if (is(se, "SummarizedExperiment0")) {
            se <- as(se, "SummarizedExperiment")
        } else if (is(se, "SummarizedExperiment")) {
            # Only to help transition from SummarizedExperiment to new
            # RangedSummarizedExperiment objects, remove once transition is
            # complete.
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
    if (run$analysis == "srnaseq") {
        bcb <- read_smallrna_counts(bcb)
    }
    bcb
}



#' Load bcbio-nextgen run
#'
#' Simply point to the final upload directory output by
#' [bcbio-nextgen](https://bcbio-nextgen.readthedocs.io/), and this function
#' will take care of the rest. It automatically imports RNA-seq counts,
#' metadata, and program versions used.
#'
#' When working in RStudio, we recommend connecting to the bcbio-nextgen run
#' directory as a remote connection over
#' [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh
#' @author Lorena Patano
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param analysis Analysis type (e.g. `rnaseq` or `srnaseq`).
#' @param groups_of_interest Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param metadata_file Custom metadata file to import. Otherwise defaults to
#'   sample metadata saved in the YAML file.
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes and normally does not need to be set.
#' @param ... Additional optional parameters passed to [SummarizedExperiment].
#'
#' @return [bcbioRnaDataSet] containing counts and metadata.
#' @export
load_run <- function(
    upload_dir = "final",
    analysis = "rnaseq",
    groups_of_interest = "description",
    metadata_file = NULL,
    organism = NULL,
    ...) {
    run <- load_run_as_list(
        upload_dir,
        analysis = analysis,
        groups_of_interest = groups_of_interest,
        organism = organism)

    # colData
    message("Processing sample metadata (colData)...")
    if (!is.null(metadata_file)) {
        colData <- read_metadata(metadata_file, lanes = run$lanes)
        run$metadata_file <- metadata_file
    } else {
        colData <- read_bcbio_metadata(run)
    }
    colData <- DataFrame(colData)

    # Subset the sample_dirs by the metadata data frame
    run$sample_dirs <- run$sample_dirs[colData$description]
    message(paste(length(run$sample_dirs), "samples matched by metadata"))

    if (analysis == "rnaseq") {
        # Metrics ----
        message("Reading metrics...")
        run$metrics <- read_bcbio_metrics(run)

        # tx2gene ----
        # Use the tx2gene file output by `bcbio-nextgen`. Alternatively,
        # we can handle this directly in R using the Ensembl annotations
        # obtained with biomaRt instead, if the file doesn't exist. This
        # is currently the case with the fast RNA-seq pipeline, for example.
        message("Processing tx2gene annotations...")
        run$tx2gene <- read_bcbio_file(
            run, file = "tx2gene.csv", col_names = FALSE)
        if (is.null(run$tx2gene)) {
            run$tx2gene <- tx2gene(run)
        }

        # Ensembl annotations ----
        message("Downloading gene annotations from Ensembl (rowData)...")
        run$ensembl_version <- listMarts() %>% .[1, 2]
        ensembl <- ensembl(run) %>% gene_level_annotations
    }

    # Read counts ----
    # Lightweight-aligned counts with tximport
    if (analysis %in% c("rnaseq", "srnaseq")) {
        # [fix] Method for loading aligned accounts (e.g. STAR)? If we add
        # support for this, we don't want to slot into txi.
        run$txi <- read_bcbio_counts(run)
    } else {
        stop("Unsupported analysis method")
    }

    # rowData ----
    # Here we need to subset the Ensembl annotations to match the genes in the
    # count matrix data
    run$txi$counts %>% nrow
    ensembl %>% nrow
    # This fails...
    all(rownames(run$txi$counts) %in% rownames(ensembl))
    # How to expand to match??

    se <- SummarizedExperiment(
        assays = SimpleList(counts = run$txi$counts),
        colData = colData,
        # rowData = rowData,
        metadata = run,
        ...)
    .bcbioRnaDataSet(se, run)
}
