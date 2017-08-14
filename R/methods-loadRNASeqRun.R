#' Load bcbio RNA-Seq Run
#'
#' Simply point to the final upload directory output by
#' [bcbio](https://bcbio-nextgen.readthedocs.io/), and this function will take
#' care of the rest. It automatically imports RNA-seq counts, metadata, and
#' program versions used.
#'
#' @rdname loadRNASeqRun
#' @name loadRNASeqRun
#'
#' @param object Path to final upload directory. This path is set when running
#'   `bcbio_nextgen -w template`.
#' @param interestingGroups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param sampleMetadataFile *Optional*. Custom metadata file containing
#'   sample information. Otherwise defaults to sample metadata saved in the YAML
#'   file.
#' @param maxSamples Maximum number of samples to calculate [DESeq2::rlog()] and
#'   [DESeq2::varianceStabilizingTransformation()] matrix. See Details.
#' @param ... Additional arguments, slotted into the [metadata()] accessor.
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @details
#' When number of samples is bigger than `maxSamples`, `rlog` and `vst` slot
#' in [SummarizedExperiment::SummarizedExperiment] will be the output of
#' [edgeR] normalization method.
#'
#' @return [bcbioRNADataSet].
#'
#' @examples
#' extraDir <- system.file("extra", package = "bcbioRnaseq")
#' uploadDir <- file.path(extraDir, "bcbio")
#' sampleMetadataFile <- file.path(extraDir, "sample_metadata.csv")
#' bcb <- loadRNASeqRun(
#'     uploadDir,
#'     interestingGroups = "group",
#'     sampleMetadataFile = sampleMetadataFile)
#' # bcbSmall <- bcb[1:999, 1:3]
NULL



# Methods ====
#' @rdname loadRNASeqRun
#' @export
setMethod("loadRNASeqRun", "character", function(
    object,
    interestingGroups = "sampleName",
    sampleMetadataFile = NULL,
    maxSamples = 50,
    ...) {
    uploadDir <- object

    # Directory paths ====
    # Check connection to final upload directory
    if (!dir.exists(uploadDir)) {
        stop("Final upload directory failed to load")
    }
    uploadDir <- normalizePath(uploadDir)
    # Find most recent nested projectDir (normally only 1)
    projectDir <- dir(uploadDir,
                      pattern = projectDirPattern,
                      full.names = FALSE,
                      recursive = FALSE)
    if (length(projectDir) != 1L) {
        stop("Uncertain about project directory location")
    }
    message(projectDir)
    match <- str_match(projectDir, projectDirPattern)
    runDate <- match[[2L]] %>% as.Date
    template <- match[[3L]]
    projectDir <- file.path(uploadDir, projectDir)

    # Project summary YAML ====
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    if (!file.exists(yamlFile)) {
        stop("YAML project summary missing")
    }
    yaml <- readYAML(yamlFile)

    # Sample directories ====
    sampleDirs <- .sampleDirs(uploadDir)

    # Sequencing lanes ====
    lanePattern <- "_L(\\d{3})"
    if (any(str_detect(sampleDirs, lanePattern))) {
        lanes <- str_match(names(sampleDirs), lanePattern) %>%
            .[, 2L] %>%
            unique %>%
            length
        message(paste(
            lanes, "sequencing lane detected", "(technical replicates)"))
    } else {
        lanes <- 1L
    }

    # Sample metadata (colData) ====
    sampleMetadata <-
        .readSampleMetadataFile(sampleMetadataFile, lanes = lanes)
    if (is.null(sampleMetadata)) {
        sampleMetadata <- .sampleYAMLMetadata(yaml)
    }
    if (!all(sampleMetadata[["sampleID"]] %in% names(sampleDirs))) {
        stop("Sample name mismatch", call. = FALSE)
    }

    # Subset sample directories by metadata ====
    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (length(sampleMetadata[["sampleID"]]) < length(sampleDirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(sampleDirs) %in% sampleMetadata[["sampleID"]]]
        message(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Genome ====
    # Use the genome build of the first sample to match
    genomeBuild <- yaml[["samples"]][[1L]][["genome_build"]]
    organism <- detectOrganism(genomeBuild)
    message(paste0("Genome: ", organism, " (", genomeBuild, ")"))
    annotable <- annotable(genomeBuild)
    tx2gene <- .tx2gene(projectDir, genomeBuild)

    # Sample metrics ====
    # Note that sample metrics used for QC plots are not currently generated
    # when using fast RNA-seq workflow. This depends upon MultiQC and aligned
    # counts generated with STAR.
    metrics <- .sampleYAMLMetrics(yaml)

    # bcbio-nextgen run information ====
    message("Reading bcbio run information")
    dataVersions <- .dataVersions(projectDir)
    programs <- .programs(projectDir)
    bcbioLog <-
        .logFile(file.path(projectDir, "bcbio-nextgen.log"))
    bcbioCommandsLog <-
        .logFile(file.path(projectDir, "bcbio-nextgen-commands.log"))

    # Metadata ====
    metadata <- list(
        version = packageVersion("bcbioRnaseq"),
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        projectDir = projectDir,
        template = template,
        runDate = runDate,
        interestingGroups = interestingGroups,
        genomeBuild = genomeBuild,
        organism = organism,
        annotable = annotable,
        tx2gene = tx2gene,
        lanes = lanes,
        yamlFile = yamlFile,
        yaml = yaml,
        metrics = metrics,
        sampleMetadataFile = sampleMetadataFile,
        dataVersions = dataVersions,
        programs = programs,
        bcbioLog = bcbioLog,
        bcbioCommandsLog = bcbioCommandsLog,
        allSamples = allSamples)
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0L) {
        metadata <- c(metadata, dots)
    }

    # tximport ====
    txi <- .tximport(sampleDirs, tx2gene = tx2gene)
    rawCounts <- txi[["counts"]]
    tmm <- .tmm(rawCounts)
    tpm <- txi[["abundance"]]

    # DESeqDataSet ====
    message("Generating internal DESeqDataSet for quality control")
    dds <- DESeqDataSetFromTximport(
        txi = txi,
        colData = sampleMetadata,
        design = formula(~1L)) %>%
        DESeq
    normalizedCounts <- counts(dds, normalized = TRUE)

    # rlog & variance
    if(nrow(sampleMetadata) > maxSamples){
        message("Data to big, skipping vst/rlog")
        rlog <- DESeqTransform(
            SummarizedExperiment(assays = log2(tmm + 1L),
                                 colData = colData(dds)))
        vst <- DESeqTransform(
            SummarizedExperiment(assays = log2(tmm + 1L),
                                 colData = colData(dds)))
    }else{
        message("Performing rlog transformation")
        rlog <- rlog(dds)
        message("Performing variance stabilizing transformation")
        vst <- varianceStabilizingTransformation(dds)
    }

    # featureCounts ====
    # STAR aligned counts, used for summary metrics. Not generated by
    # fast RNA-seq workflow.
    fcFile <- file.path(projectDir, "combined.counts")
    if (file.exists(fcFile)) {
        message("Reading STAR featureCounts aligned counts")
        fc <- read_tsv(fcFile) %>%
            as.data.frame %>%
            camel %>%
            column_to_rownames("id") %>%
            as.matrix
        if (!identical(colnames(rawCounts), colnames(fc))) {
            # Look for column name mismatch and attempt fix.
            # This is an error fix for the current bcb example dataset.
            # Safe to remove in a future update.
            # Subset columns by matching STAR sample name in metrics.
            fc <- fc %>%
                .[, camel(pull(metrics, "name"))] %>%
                # Ensure column names match tximport
                set_colnames(colnames(rawCounts))
        }
    } else {
        fc <- NULL
    }

    # Package SummarizedExperiment ====
    se <- packageSE(
        SimpleList(
            raw = rawCounts,
            normalized = normalizedCounts,
            tpm = tpm,
            tmm = tmm,
            rlog = rlog,
            vst = vst),
        colData = sampleMetadata,
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
})
