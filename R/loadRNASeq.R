#' Load bcbio RNA-Seq Data
#'
#' Simply point to the final upload directory output by
#' [bcbio](https://bcbio-nextgen.readthedocs.io/), and this function will take
#' care of the rest. It automatically imports RNA-seq counts, metadata, and
#' program versions used.
#'
#' @author Michael Steinbaugh, Lorena Pantano
#'
#' @importFrom basejump annotable camel prepareAnnotable
#'   prepareSummarizedExperiment
#' @importFrom DESeq2 DESeq DESeqDataSetFromTximport DESeqTransform rlog
#'  varianceStabilizingTransformation
#' @importFrom dplyr pull
#' @importFrom magrittr set_colnames
#' @importFrom stats formula
#' @importFrom stringr str_match
#' @importFrom tibble column_to_rownames
#' @importFrom utils packageVersion
#'
#' @param uploadDir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param interestingGroups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param sampleMetadataFile *Optional*. Custom metadata file containing
#'   sample information. Otherwise defaults to sample metadata saved in the YAML
#'   file.
#' @param maxSamples Maximum number of samples to calculate [DESeq2::rlog()] and
#'   [DESeq2::varianceStabilizingTransformation()] matrix. See Details.
#' @param annotable *Optional*. User-defined gene annotations (a.k.a.
#'   "annotable"), which will be slotted into [rowData()]. Typically this should
#'   be left undefined. By default, the function will automatically generate an
#'   annotable from the annotations available on Ensembl. If set `NULL`, then
#'   [rowData()] inside the resulting [bcbioRNASeq] object will be left empty.
#'   This is recommended for projects dealing with genes or transcripts that are
#'   poorly annotated.
#' @param ensemblVersion *Optional*. Ensembl release version. If `NULL`,
#'   defaults to current release, and does not typically need to be
#'   user-defined. This parameter can be useful for matching Ensembl annotations
#'   against an outdated bcbio annotation build.
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
#' @return [bcbioRNASeq].
#' @export
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
#' bcb <- loadRNASeq(uploadDir, interestingGroups = "group")
loadRNASeq <- function(
    uploadDir,
    interestingGroups = "sampleName",
    sampleMetadataFile = NULL,
    maxSamples = 50,
    annotable,
    ensemblVersion = NULL,
    ...) {
    # Directory paths ====
    # Check connection to final upload directory
    if (!dir.exists(uploadDir)) {
        stop("Final upload directory failed to load", call. = FALSE)
    }
    uploadDir <- normalizePath(uploadDir)
    projectDir <- dir(
        uploadDir,
        pattern = projectDirPattern,
        full.names = FALSE,
        recursive = FALSE)
    if (length(projectDir) != 1) {
        stop("Uncertain about project directory location", call. = FALSE)
    }
    message(projectDir)
    match <- str_match(projectDir, projectDirPattern)
    runDate <- match[[2]] %>%
        as.Date()
    template <- match[[3]]
    projectDir <- file.path(uploadDir, projectDir)
    sampleDirs <- .sampleDirs(uploadDir)

    # Sequencing lanes ====
    lanePattern <- "_L(\\d{3})"
    if (any(grepl(x = sampleDirs, pattern = lanePattern))) {
        lanes <- str_match(names(sampleDirs), lanePattern) %>%
            .[, 2] %>%
            unique() %>%
            length()
        message(paste(
            lanes, "sequencing lane detected", "(technical replicates)"))
    } else {
        lanes <- 1
    }

    # Project summary YAML ====
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    if (!file.exists(yamlFile)) {
        stop("'project-summary.yaml' file missing", call. = FALSE)
    }
    yaml <- readYAML(yamlFile)

    # Sample metadata ====
    if (!is.null(sampleMetadataFile)) {
        sampleMetadataFile <- normalizePath(sampleMetadataFile)
        sampleMetadata <- readSampleMetadataFile(
            file = sampleMetadataFile,
            lanes = lanes)
    } else {
        sampleMetadata <- sampleYAMLMetadata(yaml)
    }
    if (!all(sampleMetadata[["sampleID"]] %in% names(sampleDirs))) {
        stop("Sample name mismatch", call. = FALSE)
    }

    # Interesting groups ====
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    # Default to `sampleName`
    if (is.null(interestingGroups)) {
        warning(paste(
            "'interestingGroups' is 'NULL'.",
            "Defaulting to 'sampleName'."
            ), call. = FALSE)
        interestingGroups <- "sampleName"
    }
    # Check to ensure interesting groups are defined
    if (!all(interestingGroups %in% colnames(sampleMetadata))) {
        stop("Interesting groups missing in sample metadata", call. = FALSE)
    }

    # Subset sample directories by metadata ====
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
    genomeBuild <- yaml[["samples"]][[1]][["genome_build"]]
    organism <- detectOrganism(genomeBuild)
    message(paste0("Genome: ", organism, " (", genomeBuild, ")"))

    # Gene and transcript annotations ====
    if (missing(annotable)) {
        annotable <- annotable(genomeBuild, release = ensemblVersion)
    } else if (!is.null(annotable)) {
        annotable <- prepareAnnotable(annotable)
    }
    tx2gene <- .tx2gene(
        projectDir,
        organism = organism,
        release = ensemblVersion)

    # Sample metrics ====
    # Note that sample metrics used for QC plots are not currently generated
    # when using fast RNA-seq workflow. This depends upon MultiQC and aligned
    # counts generated with STAR.
    metrics <- sampleYAMLMetrics(yaml)

    # bcbio-nextgen run information ====
    message("Reading bcbio run information")
    dataVersions <- readDataVersions(
        file.path(projectDir, "data_versions.csv"))
    programs <- readProgramVersions(
        file.path(projectDir, "programs.txt"))
    bcbioLog <- readLogFile(
        file.path(projectDir, "bcbio-nextgen.log"))
    bcbioCommandsLog <- readLogFile(
        file.path(projectDir, "bcbio-nextgen-commands.log"))

    # Metadata ====
    metadata <- list(
        version = packageVersion("bcbioRNASeq"),
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        projectDir = projectDir,
        template = template,
        runDate = runDate,
        interestingGroups = interestingGroups,
        organism = organism,
        genomeBuild = genomeBuild,
        ensemblVersion = ensemblVersion,
        annotable = annotable,
        tx2gene = tx2gene,
        lanes = lanes,
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
    if (length(dots) > 0) {
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
        design = formula(~1))
    dds <- suppressWarnings(DESeq(dds))
    normalizedCounts <- counts(dds, normalized = TRUE)

    # rlog & variance ====
    if (nrow(sampleMetadata) > maxSamples) {
        warning("Dataset too large, skipping count transformations",
                call. = FALSE)
        rlog <- SummarizedExperiment(
                assays = log2(tmm + 1),
                colData = colData(dds)) %>%
            DESeqTransform()
        vst <- SummarizedExperiment(
                assays = log2(tmm + 1),
                colData = colData(dds)) %>%
            DESeqTransform()
    } else {
        message("Performing rlog transformation")
        rlog <- rlog(dds)
        message("Performing variance stabilizing transformation")
        vst <- varianceStabilizingTransformation(dds)
    }

    # STAR featureCounts ====
    # Aligned counts, used for summary metrics. Not generated for fast RNA-seq.
    fcFile <- file.path(projectDir, "combined.counts")
    if (file.exists(fcFile)) {
        message("Reading STAR featureCounts aligned counts")
        fc <- read_tsv(fcFile) %>%
            as.data.frame() %>%
            # Sanitize sampleIDs in colnames into valid names
            set_colnames(
                gsub(x = make.names(colnames(.), unique = TRUE),
                     pattern = "\\.",
                     replacement = "_")
            ) %>%
            column_to_rownames("id") %>%
            as.matrix()
        if (!identical(colnames(rawCounts), colnames(fc))) {
            # Look for column name mismatch and attempt fix.
            # This is an error fix for the current bcb example dataset.
            # Safe to remove in a future update.
            # Subset columns by matching STAR sample name in metrics.
            fc <- fc %>%
                .[, gsub(x = make.names(pull(metrics, "name"), unique = TRUE),
                        pattern = "\\.",
                        replacement = "_"), drop = FALSE] %>%
                # Ensure column names match tximport
                set_colnames(colnames(rawCounts))
        }
    } else {
        fc <- NULL
    }

    # Prepare SummarizedExperiment ====
    se <- prepareSummarizedExperiment(
        assays = list(
            raw = rawCounts,
            normalized = normalizedCounts,
            tpm = tpm,
            tmm = tmm,
            rlog = rlog,
            vst = vst),
        rowData = annotable,
        colData = sampleMetadata,
        metadata = metadata)

    # bcbioRNASeq ====
    bcb <- new("bcbioRNASeq", se)
    # Slot additional data
    bcbio(bcb, "tximport") <- txi
    bcbio(bcb, "DESeqDataSet") <- dds
    # Slot STAR featureCounts matrix, if present
    if (is.matrix(fc)) {
        bcbio(bcb, "featureCounts") <- fc
    }
    bcb
}
