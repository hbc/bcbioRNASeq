#' Load bcbio RNA-Seq Data
#'
#' Simply point to the final upload directory output by
#' [bcbio](https://bcbio-nextgen.readthedocs.io/), and this function will take
#' care of the rest. It automatically imports RNA-seq counts, metadata, and the
#' program versions used.
#'
#' When the number of samples is bigger than the `transformationLimit`, `rlog`
#' and `vst` counts will not be slotted into `assays()`. In this case, we
#' recommend visualization using `tmm` counts, which are automatically
#' calculated by [edgeR::tmm()].
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Lorena Pantano
#'
#' @importFrom basejump ensembl readYAML transcripts
#' @importFrom bcbioBase prepareSummarizedExperiment readDataVersions
#'   readLogFile readProgramVersions readSampleMetadataFile sampleDirs
#'   sampleYAMLMetadata sampleYAMLMetrics
#' @importFrom DESeq2 DESeq DESeqDataSetFromTximport DESeqTransform rlog
#'  varianceStabilizingTransformation
#' @importFrom dplyr mutate_all pull
#' @importFrom stats formula
#' @importFrom stringr str_match str_trunc
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom utils packageVersion
#'
#' @param uploadDir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param level Import counts as "`genes`" (default) or "`transcripts`".
#' @param caller Expression caller. Supports "`salmon`" (default), "`kallisto`",
#'   or "`sailfish`".
#' @param interestingGroups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param samples *Optional.* Specify a subset of samples to load. The names
#'   must match the `description` specified in the bcbio YAML metadata. If a
#'   `sampleMetadataFile` is provided, that will take priority for sample
#'   selection. Typically this can be left `NULL`.
#' @param sampleMetadataFile *Optional.* Custom metadata file containing
#'   sample information. Otherwise defaults to sample metadata saved in the YAML
#'   file. Remote URLs are supported. Typically this can be left `NULL`.
#' @param rowRanges *Required.* Genomic ranges (`GRanges`) that match the
#'   rownames of the counts matrix. If left missing (default), the annotations
#'   will be obtained automatically fron Ensembl using AnnotationHub and
#'   ensembldb. These values are also accessible with the [rowData()] function.
#' @param isSpike Genes or transcripts corresponding to FASTA spike-in
#'   sequences (e.g. ERCCs, EGFP, TDTOMATO).
#' @param organism *Optional.* Organism name. Use the full latin name (e.g.
#'   "Homo sapiens"). Normally this can be left `NULL`, and the function will
#'   detect the organism from the gene identifiers present in the rownames of
#'   the counts matrix using [detectOrganism()]. Alternatively, automatic
#'   organism detection can be overriden here, and the argument will be passed
#'   to AnnotationHub.
#' @param ensemblRelease *Optional.* Ensembl release version. If `NULL`,
#'   defaults to current release, and does not typically need to be
#'   user-defined. This parameter can be useful for matching Ensembl annotations
#'   against an outdated bcbio annotation build.
#' @param genomeBuild *Optional.* Genome build. Normally this can be left `NULL`
#'   and the build will be detected from the bcbio run data. This can be set
#'   manually (e.g. "GRCh37" for the older *Homo sapiens* reference genome).
#'   Note that this must match the genome build identifier on Ensembl for
#'   annotations to download correctly.
#' @param design DESeq2 design formula. Empty by default. Can be updated after
#'   initial data loading using the [design()] function.
#' @param transformationLimit Maximum number of samples to calculate
#'   [DESeq2::rlog()] and [DESeq2::varianceStabilizingTransformation()] matrix.
#'   It is not generally recommended to change this value. For large datasets,
#'   DESeq2 will take a really long time applying variance stabilization. See
#'   Details. Use `Inf` to always apply transformations and `0` to always skip.
#' @param ... Additional arguments, slotted into the [metadata()] accessor.
#'
#' @return `bcbioRNASeq`.
#' @export
#'
#' @seealso `help("bcbioRNASeq-class", "bcbioRNASeq")`.
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
#'
#' # Gene level
#' loadRNASeq(uploadDir, level = "genes")
#'
#' # Transcript level
#' loadRNASeq(uploadDir, level = "transcripts")
loadRNASeq <- function(
    uploadDir,
    level = c("genes", "transcripts"),
    caller = c("salmon", "kallisto", "sailfish"),
    interestingGroups = "sampleName",
    samples = NULL,
    sampleMetadataFile = NULL,
    rowRanges,
    isSpike = NULL,
    organism = NULL,
    ensemblRelease = NULL,
    genomeBuild = NULL,
    design = formula(~1),
    transformationLimit = 50L,
    ...
) {
    dots <- list(...)

    # Legacy argument matches ==================================================
    call <- match.call(expand.dots = TRUE)
    # rowRanges : annotable
    if ("annotable" %in% names(call)) {
        warn("Use `rowRanges` instead of `annotable`")
        rowRanges <- call[["annotable"]]
        dots[["annotable"]] <- NULL
    }
    # ensemblRelease : ensemblVersion
    if ("ensemblVersion" %in% names(call)) {
        warn("Use `ensemblRelease` instead of `ensemblVersion`")
        ensemblRelease <- call[["ensemblVersion"]]
        dots[["ensemblVersion"]] <- NULL
    }
    dots <- Filter(Negate(is.null), dots)

    # Assert checks ============================================================
    assert_all_are_dirs(uploadDir)
    level <- match.arg(level)
    caller <- match.arg(caller)
    assert_is_character(interestingGroups)
    assertIsAStringOrNULL(sampleMetadataFile)
    assertIsCharacterOrNULL(samples)
    if (!missing(rowRanges)) {
        assert_is_all_of(rowRanges, "GRanges")
    }
    assertIsCharacterOrNULL(isSpike)
    assertIsAStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsAStringOrNULL(genomeBuild)
    assert_is_a_number(transformationLimit)
    assert_all_are_non_negative(transformationLimit)
    assert_is_formula(design)

    # Directory paths ==========================================================
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    projectDir <- list.files(
        path = uploadDir,
        pattern = bcbioBase::projectDirPattern,
        full.names = FALSE,
        recursive = FALSE
    )
    assert_is_a_string(projectDir)
    inform(projectDir)
    match <- str_match(projectDir, bcbioBase::projectDirPattern)
    runDate <- as.Date(match[[2L]])
    template <- match[[3L]]
    projectDir <- file.path(uploadDir, projectDir)
    assert_all_are_dirs(projectDir)
    sampleDirs <- sampleDirs(uploadDir)
    assert_all_are_dirs(sampleDirs)

    # Sequencing lanes =========================================================
    if (any(grepl(x = sampleDirs, pattern = lanePattern))) {
        lanes <- str_match(names(sampleDirs), lanePattern) %>%
            .[, 2L] %>%
            unique() %>%
            length()
        inform(paste(
            lanes, "sequencing lane detected", "(technical replicates)"
        ))
    } else {
        lanes <- 1L
    }
    assert_is_an_integer(lanes)

    # Project summary YAML =====================================================
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    assert_all_are_existing_files(yamlFile)
    yaml <- readYAML(yamlFile)
    assert_is_list(yaml)

    # Column data ==============================================================
    if (is_a_string(sampleMetadataFile)) {
        colData <- readSampleMetadataFile(sampleMetadataFile, lanes = lanes)
    } else {
        colData <- sampleYAMLMetadata(yaml)
        if (is.character(samples)) {
            assert_is_subset(samples, colData[["description"]])
            colData <- colData %>%
                .[which(samples %in% .[["description"]]), , drop = FALSE]
        }
    }

    # Interesting groups =======================================================
    interestingGroups <- camel(interestingGroups)
    assert_is_subset(interestingGroups, colnames(colData))

    # Subset sample directories by metadata ====================================
    samples <- colData[["sampleID"]]
    assert_are_intersecting_sets(samples, names(sampleDirs))
    if (length(samples) < length(sampleDirs)) {
        inform(paste(
            "Loading a subset of samples:",
            str_trunc(toString(samples), width = 80L),
            sep = "\n"
        ))
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(.) %in% samples]
    } else {
        allSamples <- TRUE
    }

    # Genome ===================================================================
    # Genome build
    if (!is_a_string(genomeBuild)) {
        # Detect from the bcbio project summary YAML
        genomeBuild <- yaml %>%
            .[["samples"]] %>%
            .[[1L]] %>%
            .[["genome_build"]]
    }
    assert_is_a_string(genomeBuild)

    # Organism
    if (is.null(organism) && is_a_string(genomeBuild)) {
        inform("Detecting organism from genome build")
        organism <- detectOrganism(genomeBuild)
    }
    assert_is_a_string(organism)
    inform(paste("Genome:", organism, paste0("(", genomeBuild, ")")))

    # Row data =================================================================
    if (missing(rowRanges)) {
        # ah = AnnotationHub
        ah <- ensembl(
            organism = organism,
            format = level,
            genomeBuild = genomeBuild,
            release = ensemblRelease,
            return = "GRanges",
            metadata = TRUE
        )
        assert_is_list(ah)
        assert_are_identical(names(ah), c("data", "metadata"))
        rowRanges <- ah[["data"]]
        assert_is_all_of(rowRanges, "GRanges")
        ahMeta <- ah[["metadata"]]
        assert_all_are_matching_regex(
            x = ahMeta[["id"]],
            pattern = "^AH\\d+$"
        )
    } else {
        ahMeta <- NULL
    }

    # Sample metrics ===========================================================
    # Note that sample metrics used for QC plots are not currently generated
    # when using fast RNA-seq workflow. This depends upon MultiQC and aligned
    # counts generated with STAR.
    inform("Reading sample metrics")
    metrics <- sampleYAMLMetrics(yaml)
    assert_is_data.frame(metrics)

    # bcbio-nextgen run information ============================================
    inform("Reading bcbio run information")
    dataVersions <- readDataVersions(
        file = file.path(projectDir, "data_versions.csv")
    )
    assert_is_tbl_df(dataVersions)

    programVersions <- readProgramVersions(
        file = file.path(projectDir, "programs.txt")
    )
    assert_is_tbl_df(programVersions)

    bcbioLogFile <- list.files(
        path = projectDir,
        pattern = "bcbio-nextgen.log",
        full.names = TRUE,
        recursive = FALSE
    )
    assert_is_a_string(bcbioLogFile)
    inform(basename(bcbioLogFile))
    bcbioLog <- readLogFile(bcbioLogFile)
    assert_is_character(bcbioLog)

    bcbioCommandsLogFile <- list.files(
        path = projectDir,
        pattern = "bcbio-nextgen-commands.log",
        full.names = TRUE,
        recursive = FALSE
    )
    assert_is_a_string(bcbioCommandsLogFile)
    inform(basename(bcbioCommandsLogFile))
    bcbioCommandsLog <- readLogFile(bcbioCommandsLogFile)
    assert_is_character(bcbioCommandsLog)

    # tximport =================================================================
    if (level == "transcripts") {
        txOut = TRUE
    } else {
        txOut = FALSE
    }
    # Attempt to use `tx2gene.csv` saved in project directory
    tx2gene <- .tx2gene(
        projectDir = projectDir,
        organism = organism,
        release = ensemblRelease
    )
    txi <- .tximport(
        sampleDirs = sampleDirs,
        type = caller,
        txIn = TRUE,
        txOut = txOut,
        tx2gene = tx2gene
    )
    # abundance = transcripts per million (TPM)
    # counts = raw counts
    # length = average transcript length
    # countsFromAbundance = character describing TPM
    tpm <- txi[["abundance"]]
    assert_is_matrix(tpm)
    counts <- txi[["counts"]]
    assert_is_matrix(counts)
    length <- txi[["length"]]
    assert_is_matrix(length)
    countsFromAbundance <- txi[["countsFromAbundance"]]
    assert_is_character(countsFromAbundance)

    # colData ==================================================================
    # This step is necessary for samples have been sanitized with
    # `make.names()`, which can cause samples that are now prefixed with `X` to
    # become out of order in the rows. This will cause the
    # `DESeqDataSetFromTximport` call below to error out because of a
    # `txi$counts` colnames / colData rownames mismatch.
    colData <- colData %>%
        as.data.frame() %>%
        .[colnames(counts), , drop = FALSE] %>%
        rownames_to_column() %>%
        mutate_all(as.factor) %>%
        mutate_all(droplevels) %>%
        column_to_rownames()

    # Gene-level specific calculations =========================================
    if (level == "genes") {
        inform(paste(
            "Generating internal DESeqDataSet using DESeq2",
            packageVersion("DESeq2")
        ))
        if (!is(design, "formula")) {
            design <- formula(~1)  # nolint
        }
        dds <- DESeqDataSetFromTximport(
            txi = txi,
            colData = colData,
            design = design
        )
        # Suppressing warnings here for empty design formula (`~1`)
        dds <- suppressWarnings(DESeq(dds))
        # Variance stabilizing transformations
        if (nrow(colData) > transformationLimit) {
            warn(paste(
                "Dataset contains many samples.",
                "Skipping variance stabilizing transformations."
            ))
            rlog <- NULL
            vst <- NULL
        } else {
            inform("Performing rlog transformation")
            rlog <- rlog(dds)
            inform("Performing variance stabilizing transformation")
            vst <- varianceStabilizingTransformation(dds)
        }
    } else {
        dds <- NULL
        rlog <- NULL
        vst <- NULL
    }

    # Metadata =================================================================
    metadata <- list(
        "version" = packageVersion,
        "uploadDir" = uploadDir,
        "sampleDirs" = sampleDirs,
        "projectDir" = projectDir,
        "template" = template,
        "runDate" = runDate,
        "level" = level,
        "caller" = caller,
        "interestingGroups" = interestingGroups,
        "organism" = organism,
        "genomeBuild" = genomeBuild,
        "ensemblRelease" = as.integer(ensemblRelease),
        "annotationHub" = as.list(ahMeta),
        "tx2gene" = tx2gene,
        "countsFromAbundance" = countsFromAbundance,
        "lanes" = lanes,
        "yaml" = yaml,
        "metrics" = metrics,
        "sampleMetadataFile" = as.character(sampleMetadataFile),
        "dataVersions" = dataVersions,
        "programVersions" = programVersions,
        "bcbioLog" = bcbioLog,
        "bcbioCommandsLog" = bcbioCommandsLog,
        "allSamples" = allSamples,
        "design" = design,
        "loadRNASeq" = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots)) {
        metadata <- c(metadata, dots)
    }

    # Return =============================================
    assays = list(
        "raw" = counts,
        "tpm" = tpm,
        "length" = length,
        "dds" = dds,
        "rlog" = rlog,
        "vst" = vst
    )
    rse <- prepareSummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        isSpike = isSpike
    )
    assert_is_all_of(rse, "RangedSummarizedExperiment")
    new("bcbioRNASeq", rse)
}
