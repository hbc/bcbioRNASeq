#' Load `bcbio` Run Data
#'
#' Simply point to the final upload directory output by
#' [bcbio](https://bcbio-nextgen.readthedocs.io/), and this function will take
#' care of the rest. It automatically imports RNA-seq counts, metadata, and
#' program versions used.
#'
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param analysis Analysis type. Supports RNA-seq (`rnaseq`; **default**) or
#'   small RNA-seq (`srnaseq`).
#' @param interesting_groups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param sample_metadata_file *Optional*. Custom metadata file containing
#'   sample information. Otherwise defaults to sample metadata saved in the YAML
#'   file.
#' @param ... Additional arguments, slotted into the [metadata()] accessor.
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @return [bcbioRNADataSet].
#' @export
#'
#' @examples
#' extra_dir <- system.file("extra", package = "bcbioRnaseq")
#' upload_dir <- file.path(extra_dir, "bcbio")
#' sample_metadata_file <- file.path(extra_dir, "sample_metadata.csv")
#' bcb <- load_run(
#'     upload_dir,
#'     interesting_groups = "group",
#'     sample_metadata_file = sample_metadata_file)
load_run <- function(
    upload_dir = "final",
    analysis = "rnaseq",
    interesting_groups = "sample_name",
    sample_metadata_file = NULL,
    ...) {
    # Analysis type ====
    supported_analyses <- c("rnaseq", "srnaseq")
    if (!analysis %in% supported_analyses) {
        stop(paste("Supported analyses:", toString(supported_analyses)))
    }


    # Directory paths ====
    # Check connection to final upload directory
    if (!dir.exists(upload_dir)) {
        stop("Final upload directory failed to load")
    }
    upload_dir <- normalizePath(upload_dir)
    # Find most recent nested project_dir (normally only 1)
    project_dir <- dir(upload_dir,
                       pattern = project_dir_pattern,
                       full.names = FALSE,
                       recursive = FALSE)
    if (length(project_dir) != 1L) {
        stop("Uncertain about project directory location")
    }
    message(project_dir)
    match <- str_match(project_dir, project_dir_pattern)
    run_date <- match[[2L]] %>% as.Date
    template <- match[[3L]]
    project_dir <- file.path(upload_dir, project_dir)


    # Project summary YAML ====
    yaml_file <- file.path(project_dir, "project-summary.yaml")
    if (!file.exists(yaml_file)) {
        stop("YAML project summary missing")
    }
    yaml <- read_yaml(yaml_file)


    # Sample directories ====
    sample_dirs <- .sample_dirs(upload_dir)


    # Sequencing lanes ====
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2L] %>%
            unique %>%
            length
        message(paste(
            lanes, "sequencing lane detected", "(technical replicates)"))
    } else {
        lanes <- 1L
    }


    # Sample metadata (colData) ====
    sample_metadata <-
        .sample_metadata_file(sample_metadata_file, lanes = lanes)
    if (is.null(sample_metadata)) {
        sample_metadata <- .sample_yaml_metadata(yaml)
    }


    # Subset sample directories by metadata ====
    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (length(sample_metadata[["sample_id"]]) < length(sample_dirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        all_samples <- FALSE
        sample_dirs <- sample_dirs %>%
            .[names(sample_dirs) %in% sample_metadata[["sample_id"]]]
        message(paste(length(sample_dirs), "samples matched by metadata"))
    } else {
        all_samples <- TRUE
    }


    # Genome ====
    # Use the genome build of the first sample to match
    genome_build <- yaml[["samples"]][[1L]][["genome_build"]]
    organism <- .detect_organism(genome_build)
    message(paste("Genome:", organism, genome_build))
    annotable <- annotable(genome_build)
    tx2gene <- .tx2gene(project_dir, genome_build)


    # Sample metrics ====
    # Note that sample metrics used for QC plots are not currently generated
    # when using fast RNA-seq workflow. This depends upon MultiQC and aligned
    # counts generated with STAR.
    metrics <- .sample_yaml_metrics(yaml)


    # bcbio-nextgen run information ====
    message("Reading bcbio run information")
    data_versions <- .data_versions(project_dir)
    programs <- .programs(project_dir)
    bcbio_nextgen_log <-
        .log_file(file.path(project_dir, "bcbio-nextgen.log"))
    bcbio_nextgen_commands_log <-
        .log_file(file.path(project_dir, "bcbio-nextgen-commands.log"))


    # Metadata ====
    metadata <- list(
        analysis = analysis,
        upload_dir = upload_dir,
        sample_dirs = sample_dirs,
        project_dir = project_dir,
        template = template,
        run_date = run_date,
        interesting_groups = interesting_groups,
        organism = organism,
        genome_build = genome_build,
        ensembl_version = ensembl_version(),
        annotable = annotable,
        tx2gene = tx2gene,
        lanes = lanes,
        yaml_file = yaml_file,
        yaml = yaml,
        metrics = metrics,
        sample_metadata_file = sample_metadata_file,
        data_versions = data_versions,
        programs = programs,
        bcbio_nextgen_log = bcbio_nextgen_log,
        bcbio_nextgen_commands_log = bcbio_nextgen_commands_log,
        all_samples = all_samples)

    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0L) {
        metadata <- c(metadata, dots)
    }


    # tximport ====
    txi <- .tximport(sample_dirs, tx2gene = tx2gene)
    raw_counts <- txi[["counts"]]
    tpm <- txi[["abundance"]]


    # DESeqDataSet ====
    message("Generating internal DESeqDataSet for quality control")
    dds <- DESeqDataSetFromTximport(
        txi = txi,
        colData = sample_metadata,
        design = formula(~1L)) %>%
        DESeq


    # featureCounts ====
    # STAR aligned counts, used for summary metrics. Not generated by
    # fast RNA-seq workflow.
    fc_file <- file.path(project_dir, "combined.counts")
    if (file.exists(fc_file)) {
        message("Reading STAR/featureCounts aligned counts")
        fc <- read_tsv(fc_file) %>%
            as.data.frame %>%
            snake(rownames = FALSE) %>%
            column_to_rownames("id") %>%
            as.matrix
        if (!identical(colnames(raw_counts), colnames(fc))) {
            # Look for column name mismatch and attempt fix.
            # This is an error fix for the current bcb example dataset.
            # Safe to remove in a future update.
            # Subset columns by matching STAR sample name in metrics.
            fc <- fc %>%
                .[, snake(pull(metrics, "name"))] %>%
                # Ensure column names match tximport
                set_colnames(pull(metrics, "sample_name"))
        }
    } else {
        fc <- NULL
    }


    # Package SummarizedExperiment ====
    se <- packageSE(
        assays = SimpleList(
            raw_counts = raw_counts,
            tpm = tpm,
            normalized_counts = counts(dds, normalized = TRUE),
            tmm = .tmm(raw_counts),
            rlog = rlog(dds),
            vst = varianceStabilizingTransformation(dds)),
        colData = sample_metadata,
        rowData = annotable,
        metadata = metadata)


    # bcbioRNADataSet ====
    bcb <- new("bcbioRNADataSet", se)
    # Slot additional callers
    bcbio(bcb, "tximport") <- txi
    bcbio(bcb, "DESeqDataSet") <- dds
    if (is.matrix(fc)) {
        bcbio(bcb, "featureCounts") <- fc
    }
    bcb
}
