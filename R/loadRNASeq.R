#' Load bcbio RNA-Seq Data
#'
#' Simply point to the final upload directory output by
#' [bcbio](https://bcbio-nextgen.readthedocs.io/), and this function will take
#' care of the rest. It automatically imports RNA-seq counts, metadata, and the
#' program versions used.
#'
#' When the number of samples is bigger than the `transformationLimit`, `rlog`
#' and `vst` counts will not be slotted into `assays()`. In this case, we
#' recommend visualization using [tmm()] counts, which are automatically
#' calculated using edgeR.
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Lorena Pantano
#'
#' @inheritParams general
#' @param uploadDir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param level Import counts as "`genes`" (default) or "`transcripts`".
#' @param caller Expression caller. Supports "`salmon`" (default), "`kallisto`",
#'   or "`sailfish`".
#' @param samples *Optional.* Specify a subset of samples to load. The names
#'   must match the `description` specified in the bcbio YAML metadata. If a
#'   `sampleMetadataFile` is provided, that will take priority for sample
#'   selection. Typically this can be left unset.
#' @param sampleMetadataFile *Optional.* Custom metadata file containing
#'   sample information. Otherwise defaults to sample metadata saved in the YAML
#'   file. Remote URLs are supported. Typically this can be left unset.
#' @param organism Organism name. Use the full latin name (e.g.
#'   "Homo sapiens"), since this will be input downstream to
#'   AnnotationHub and ensembldb, unless `gffFile` is set.
#' @param genomeBuild *Optional.* Ensembl genome build name (e.g. "GRCh38").
#'   This will be passed to AnnotationHub for `EnsDb` annotation matching,
#'   unless `gffFile` is set.
#' @param ensemblRelease *Optional.* Ensembl release version. If unset,
#'   defaults to current release, and does not typically need to be
#'   user-defined. Passed to AnnotationHub for `EnsDb` annotation matching,
#'   unless `gffFile` is set.
#' @param isSpike *Optional.* Gene names corresponding to FASTA spike-in
#'   sequences (e.g. ERCCs, EGFP, TDTOMATO).
#' @param gffFile *Optional, not recommended.* By default, we recommend leaving
#'   this `NULL` for genomes that are supported on Ensembl. In this case, the
#'   row annotations ([rowRanges()]) will be obtained automatically from Ensembl
#'   by passing the `organism`, `genomeBuild`, and `ensemblRelease` arguments to
#'   AnnotationHub and ensembldb. For a genome that is not supported on Ensembl
#'   and/or AnnotationHub, a GFF/GTF (General Feature Format) file is required.
#'   Generally, we recommend using a GTF (GFFv2) file here over a GFF3 file if
#'   possible, although all GFF formats are supported. The function will
#'   internally generate a `TxDb` containing transcript-to-gene mappings and
#'   construct a `GRanges` object containing the genomic ranges ([rowRanges()]).
#' @param transformationLimit Maximum number of samples to calculate
#'   [DESeq2::rlog()] and [DESeq2::varianceStabilizingTransformation()] matrix.
#'   For large datasets, DESeq2 is slow to apply variance stabilization. In this
#'   case, we recommend loading up the dataset in a high-performance computing
#'   environment. Use `Inf` to always apply and `-Inf` to always skip.
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
#' loadRNASeq(
#'     uploadDir = uploadDir,
#'     level = "genes",
#'     organism = "Mus musculus"
#' )
#'
#' # Transcript level
#' loadRNASeq(
#'     uploadDir = uploadDir,
#'     level = "transcripts",
#'     organism = "Mus musculus"
#' )
loadRNASeq <- function(
    uploadDir,
    level = c("genes", "transcripts"),
    caller = c("salmon", "kallisto", "sailfish"),
    samples = NULL,
    sampleMetadataFile = NULL,
    interestingGroups = "sampleName",
    organism,
    ensemblRelease = NULL,
    genomeBuild = NULL,
    isSpike = NULL,
    gffFile = NULL,
    transformationLimit = 50L,
    ...
) {
    dots <- list(...)

    # Legacy arguments =========================================================
    call <- match.call(expand.dots = TRUE)
    # annotable
    if ("annotable" %in% names(call)) {
        abort("`annotable` is defunct. Consider using `gffFile` instead.")
    }
    # ensemblVersion
    if ("ensemblVersion" %in% names(call)) {
        warn("Use `ensemblRelease` instead of `ensemblVersion`")
        ensemblRelease <- call[["ensemblVersion"]]
        dots[["ensemblVersion"]] <- NULL
    }
    # organism missing
    if (!"organism" %in% names(call)) {
        abort("`organism` is now required")
    }
    dots <- Filter(Negate(is.null), dots)

    # Assert checks ============================================================
    assert_all_are_dirs(uploadDir)
    level <- match.arg(level)
    caller <- match.arg(caller)
    assertIsAStringOrNULL(sampleMetadataFile)
    assertIsCharacterOrNULL(samples)
    assert_is_character(interestingGroups)
    assert_is_a_string(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsAStringOrNULL(genomeBuild)
    assertIsCharacterOrNULL(isSpike)
    assertIsAStringOrNULL(gffFile)
    if (is_a_string(gffFile)) {
        assert_all_are_existing_files(gffFile)
    }
    assert_is_a_number(transformationLimit)

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
    # Ensure all columns are factor
    colData <- sanitizeSampleData(colData)

    # Interesting groups =======================================================
    interestingGroups <- camel(interestingGroups)
    assert_is_subset(interestingGroups, colnames(colData))

    # Subset sample directories by metadata ====================================
    samples <- colData[["sampleID"]]
    assert_is_subset(samples, names(sampleDirs))
    if (length(samples) < length(sampleDirs)) {
        inform(paste(
            "Loading a subset of samples:",
            str_trunc(toString(samples), width = 80L),
            sep = "\n"
        ))
        allSamples <- FALSE
        sampleDirs <- sampleDirs[samples]
    } else {
        allSamples <- TRUE
    }

    # Row data =================================================================
    rowRangesMetadata <- NULL
    txdb <- NULL
    tx2gene <- NULL
    if (is_a_string(gffFile)) {
        txdb <- makeTxDbFromGFF(gffFile)
        rowRanges <- genes(txdb)
        # FIXME Need to sanitize rowRanges here
        # Transcript-to-gene mappings
        if (level == "transcripts") {
            transcripts <- transcripts(txdb, columns = c("tx_name", "gene_id"))
            tx2gene <- mcols(transcripts) %>%
                as.data.frame() %>%
                set_colnames(c("txID", "geneID")) %>%
                set_rownames(.[["txID"]])
        }
    } else {
        # ah: AnnotationHub
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
        rowRangesMetadata <- ah[["metadata"]]
        assert_is_data.frame(rowRangesMetadata)
    }

    # Sample metrics ===========================================================
    # Note that sample metrics used for QC plots are not currently generated
    # when using fast RNA-seq workflow. This depends upon MultiQC and aligned
    # counts generated with STAR.
    inform("Reading sample metrics")
    metrics <- sampleYAMLMetrics(yaml)
    assert_is_data.frame(metrics)

    # bcbio run information ====================================================
    dataVersions <- readDataVersions(
        file = file.path(projectDir, "data_versions.csv")
    )
    assert_is_tbl_df(dataVersions)

    programVersions <- readProgramVersions(
        file = file.path(projectDir, "programs.txt")
    )
    assert_is_tbl_df(programVersions)

    bcbioLog <- readLogFile(
        file = file.path(projectDir, "bcbio-nextgen.log")
    )
    assert_is_character(bcbioLog)

    bcbioCommandsLog <- readLogFile(
        file = file.path(projectDir, "bcbio-nextgen-commands.log")
    )
    assert_is_character(bcbioCommandsLog)

    # tximport =================================================================
    # Run this step only after the required metadata has imported successfully
    if (level == "transcripts") {
        txOut <- TRUE
    } else {
        txOut <- FALSE
    }
    # Attempt to use `tx2gene.csv` saved in project directory
    tx2gene <- .tx2gene(
        projectDir = projectDir,
        organism = organism,
        release = ensemblRelease,
        genomeBuild = genomeBuild
    )
    txi <- .tximport(
        sampleDirs = sampleDirs,
        type = caller,
        txIn = TRUE,
        txOut = txOut,
        tx2gene = tx2gene
    )
    # abundance = transcripts per million (TPM)
    # raw = raw counts
    # length = average transcript length
    # countsFromAbundance = character describing TPM
    tpm <- txi[["abundance"]]
    raw <- txi[["counts"]]
    length <- txi[["length"]]
    countsFromAbundance <- txi[["countsFromAbundance"]]
    assert_is_matrix(tpm)
    assert_is_matrix(raw)
    assert_is_matrix(length)
    assert_is_character(countsFromAbundance)

    # Ensure `colData` matches the colnames in `assays()`
    colData <- colData[colnames(raw), , drop = FALSE]

    # Gene-level variance stabilization ========================================
    rlog <- NULL
    vst <- NULL
    if (level == "genes") {
        if (nrow(colData) <= transformationLimit) {
            inform(paste(
                "Calculating variance stabilizations using DESeq2",
                packageVersion("DESeq2")
            ))
            dds <- DESeqDataSetFromTximport(
                txi = txi,
                colData = colData,
                # Use an empty design formula
                design = ~ 1  # nolint
            )
            # Suppress warning about empty design formula
            dds <- suppressWarnings(DESeq(dds))
            inform("Applying rlog transformation")
            rlog <- assay(rlog(dds))
            inform("Applying variance stabilizing transformation")
            vst <- assay(varianceStabilizingTransformation(dds))
        } else {
            warn(paste(
                "Dataset contains many samples.",
                "Skipping variance stabilizing transformations."
            ))
        }
    }

    # Assays ===================================================================
    assays <- list(
        "raw" = raw,
        "tpm" = tpm,
        "length" = length,
        "rlog" = rlog,
        "vst" = vst
    )

    # Metadata =================================================================
    metadata <- list(
        "version" = packageVersion,
        "level" = level,
        "caller" = caller,
        "countsFromAbundance" = countsFromAbundance,
        "uploadDir" = uploadDir,
        "sampleDirs" = sampleDirs,
        "sampleMetadataFile" = as.character(sampleMetadataFile),
        "projectDir" = projectDir,
        "template" = template,
        "runDate" = runDate,
        "interestingGroups" = interestingGroups,
        "organism" = organism,
        "genomeBuild" = as.character(genomeBuild),
        "ensemblRelease" = as.integer(ensemblRelease),
        "rowRangesMetadata" = rowRangesMetadata,
        "gffFile" = as.character(gffFile),
        "txdb" = txdb,
        "tx2gene" = tx2gene,
        "lanes" = lanes,
        "yaml" = yaml,
        "metrics" = metrics,
        "dataVersions" = dataVersions,
        "programVersions" = programVersions,
        "bcbioLog" = bcbioLog,
        "bcbioCommandsLog" = bcbioCommandsLog,
        "allSamples" = allSamples,
        "loadRNASeq" = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots)) {
        assert_are_disjoint_sets(metadata, dots)
        metadata <- c(metadata, dots)
    }

    # Return =============================================
    .new.bcbioRNASeq(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        isSpike = isSpike
    )
}
