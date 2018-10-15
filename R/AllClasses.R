setClassUnion(name = "missingOrNULL", members = c("missing", "NULL"))



# bcbioRNASeq ==================================================================
# TODO Improve the object documentation here.
# TODO Consider adding countsFromAbundance check against tximport formals.
#' bcbio RNA-Seq Data Set
#'
#' `bcbioRNASeq` is an S4 class that extends `RangedSummarizedExperiment`, and
#' is designed to store a [bcbio](https://bcbio-nextgen.readthedocs.org) RNA-seq
#' analysis.
#'
#' @note `bcbioRNASeq` extended `SummarizedExperiment` prior to v0.2.0, where we
#'   migrated to `RangedSummarizedExperiment`.
#'
#' @section Automatic metadata:
#'
#' The [metadata()] slot always contains:
#'
#' - Object version.
#' - bcbio data provenance information.
#' - File paths and timestamps.
#' - R session information.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh, Lorena Pantano
#' @export
#'
#' @seealso [bcbioRNASeq()].
setClass(Class = "bcbioRNASeq", contains = "RangedSummarizedExperiment")

setValidity(
    Class = "bcbioRNASeq",
    method = function(object) {
        # Return invalid for all objects older than v0.2.
        stopifnot(metadata(object)[["version"]] >= 0.2)

        assert_is_all_of(object, "RangedSummarizedExperiment")
        assert_has_dimnames(object)

        # Metadata -------------------------------------------------------------
        # Require the user to update to Bioconductor version.
        metadata <- metadata(object)

        # Detect legacy metrics.
        if (is.data.frame(metadata[["metrics"]])) {
            stop(paste(
                "`metrics` saved in `metadata()` instead of `colData()`.",
                updateMessage
            ))
        }

        # Check that interesting groups defined in metadata are valid.
        assert_is_subset(
            x = metadata[["interestingGroups"]],
            y = colnames(colData(object))
        )

        # Error on legacy slot detection.
        intersect <- intersect(
            x = names(metadata),
            y = c(
                "design",
                "devtoolsSessionInfo",
                "ensemblVersion",
                "gtf",
                "gtfFile",
                "missingGenes",
                "programs",
                "rowRangesMetadata",
                "template",
                "utilsSessionInfo",
                "yamlFile"
            )
        )
        if (has_length(intersect)) {
            stop(paste(
                paste(
                    "Legacy metadata slots:",
                    toString(sort(intersect))
                ),
                updateMessage,
                sep = "\n"
            ))
        }

        # Class checks (order independent).
        checkClasses(
            object = metadata,
            expected = list(
                allSamples = "logical",
                bcbioCommandsLog = "character",
                bcbioLog = "character",
                call = "call",
                caller = "character",
                countsFromAbundance = "character",
                dataVersions = "DataFrame",
                date = "Date",
                ensemblRelease = "integer",
                genomeBuild = "character",
                gffFile = "character",
                interestingGroups = "character",
                lanes = "integer",
                level = "character",
                organism = "character",
                programVersions = "DataFrame",
                projectDir = "character",
                runDate = "Date",
                sampleDirs = "character",
                sampleMetadataFile = "character",
                tx2gene = "Tx2Gene",
                sessionInfo = "session_info",
                uploadDir = "character",
                version = "package_version",
                wd = "character",
                yaml = "list"
            ),
            subset = TRUE
        )

        # Additional assert checks.
        assert_is_subset(metadata[["caller"]], validCallers)
        assert_is_subset(metadata[["level"]], validLevels)
        if (is.character(metadata[["countsFromAbundance"]])) {
            assert_is_subset(
                x = metadata[["countsFromAbundance"]],
                y = eval(formals(tximport)[["countsFromAbundance"]])
            )
        }

        # Assays ---------------------------------------------------------------
        assayNames <- assayNames(object)
        assert_is_subset(requiredAssays, assayNames)
        # Check that all assays are matrices.
        # Note that in previous versions, we slotted `DESeqDataSet` and
        # `DESeqTransform`, which can result in metadata mismatches because
        # those objects contain their own `colData()` and `rowData()`.
        valid <- vapply(
            X = assays(object),
            FUN = is.matrix,
            FUN.VALUE = logical(1L),
            USE.NAMES = TRUE
        )
        if (!all(valid)) {
            stop(paste(
                paste(
                    "Assays that are not matrix:",
                    toString(names(valid[!valid]))
                ),
                updateMessage,
                sep = "\n"
            ))
        }

        # Caller-specific checks.
        caller <- metadata[["caller"]]
        if (caller %in% tximportCallers) {
            assert_is_subset(tximportAssays, assayNames)
        } else if (caller %in% featureCountsCallers) {
            assert_is_subset(featureCountsAssays, assayNames)
        }

        # Check for average transcript length matrix, if necessary.
        if (
            metadata[["caller"]] %in% tximportCallers &&
            metadata[["countsFromAbundance"]] == "no"
        ) {
            assert_is_subset("avgTxLength", assayNames)
        }

        # Row data -------------------------------------------------------------
        assert_is_all_of(rowRanges(object), "GRanges")
        rowData <- rowData(object)
        if (has_length(colnames(rowData))) {
            # FIXME Consider requiring description here.
            checkClasses(
                object = rowData,
                expected = list(
                    broadClass = "factor",
                    geneBiotype = "factor",
                    # Factor needed here for transcript-level data.
                    geneID = c("character", "factor"),
                    geneName = c("character", "factor")
                ),
                subset = TRUE
            )
        }

        # Column data ----------------------------------------------------------
        colData <- colData(object)
        assert_is_subset("sampleName", colnames(colData))
        # Error on legacy metrics columns.
        assert_are_disjoint_sets(colnames(colData), legacyMetricsCols)

        TRUE
    }
)



# DESeqAnalysis ================================================================
#' DESeq2 Analysis Container
#'
#' Class containing all elements generated during differential expression
#' analysis with DESeq2. This class is essentially a `list` with validity checks
#' to ensure `DESeqTransform` and `DESeqResults` correspond to the
#' `DESeqDataSet`.
#'
#' @section DESeqDataSet:
#'
#' We recommend generating the `DESeqDataSet` by coercion from `bcbioRNASeq`
#' object using `as(dds, "bcbioRNASeq")`. Don't use the [DESeq2::DESeqDataSet()]
#' or [DESeq2::DESeqDataSetFromMatrix()] constructors to generate the
#' `DESeqDataSet` object.
#'
#' @section DESeqResults:
#'
#' Don't modify any of the `DESeqResults` objects manually. This includes
#' rearranging the rows or dropping genes without adjusted P values. We'll take
#' care of this automatically in supported methods.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @slot data `DESeqDataSet`.
#' @slot transform `DESeqTransform`.
#' @slot results `list`. One or more unshrunken `DESeqResults`.
#' @slot lfcShrink `list`. One or more shrunken `DESeqResults`.
#'
#' @seealso [DESeqAnalysis()].
setClass(
    Class = "DESeqAnalysis",
    slots = c(
        data = "DESeqDataSet",
        transform = "DESeqTransform",
        results = "list",
        lfcShrink = "list"
    ),
    prototype = list(
        lfcShrink = list()
    )
)

setValidity(
    Class = "DESeqAnalysis",
    method = function(object) {
        data <- slot(object, "data")
        transform <- slot(object, "transform")
        results <- slot(object, "results")
        lfcShrink <- slot(object, "lfcShrink")

        stopifnot(is(data, "DESeqDataSet"))
        assertHasValidDimnames(data)

        # Require gene-to-symbol mappings.
        assert_is_subset(
            x = c("geneID", "geneName"),
            y = colnames(rowData(data))
        )

        # DESeqDataSet and DESeqTransform must match.
        assert_are_identical(
            x = dimnames(data),
            y = dimnames(transform)
        )

        # DESeqDataSet and DESeqResults must match.
        invisible(lapply(
            X = results,
            FUN = function(x) {
                assert_are_identical(
                    x = rownames(x),
                    y = rownames(data)
                )
            }
        ))

        # DESeqResults and lfcShrink rownames must match, if shrinkage has been
        # calculated.
        if (length(lfcShrink) > 0L) {
            invisible(mapply(
                unshrunken = results,
                shrunken = lfcShrink,
                FUN = function(unshrunken, shrunken) {
                    assert_are_identical(
                        x = rownames(unshrunken),
                        y = rownames(shrunken)
                    )
                }
            ))
        }

        TRUE
    }
)



# DESeqResultsTables ===========================================================
#' DESeq2 Differential Expression Results Tables
#'
#' `DESeqResults` object with `DataFrame` subsets and file path information.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @slot results `DESeqResults`. Original unmodified `DESeqResults`. Should
#'   contain all genes, including those with `NA` adjusted *P* values.
#' @slot deg `list`. Differentially expressed genes. Contains `character`
#'   vectors of genes that are upregulated (`up`) or downregulated (`down`).
#'   Values map to the `rownames` of the internal `DESeqResults`. These are
#'   genes that pass `alpha` and `lfcThreshold` cutoffs set in
#'   [DESeq2::results()] call.
#' @slot counts `matrix`. Normalized counts matrix.
#' @slot rowData `DataFrame`. Row annotations.
#' @slot sampleNames `character`. Human-friendly sample names. Must contain
#'   [names()] that map to the [colnames()] of the `DESeqDataSet`.
#' @slot metadata `list`. Metadata. Contains file paths and information on
#'   whether we're writing locally or to Dropbox.
#'
#' @seealso [DESeqResultsTables()].
setClass(
    Class = "DESeqResultsTables",
    slots = c(
        results = "DESeqResults",
        deg = "list",
        counts = "matrix",
        rowData = "DataFrame",
        sampleNames = "character",
        metadata = "list"
    ),
    prototype = list(
        counts = matrix(),
        rowData = DataFrame(),
        sampleNames = character(),
        metadata = list()
    )
)

setValidity(
    Class = "DESeqResultsTables",
    method = function(object) {
        results <- slot(object, "results")
        stopifnot(is(results, "DESeqResults"))
        assert_is_a_string(contrastName(results))

        deg <- slot(object, "deg")
        up <- deg[["up"]]
        assert_is_character(up)
        assert_is_subset(up, rownames(results))
        down <- deg[["down"]]
        assert_is_character(down)
        assert_is_subset(down, rownames(results))
        assert_are_disjoint_sets(up, down)

        alpha <- metadata(results)[["alpha"]]
        assert_is_a_number(alpha)
        lfcThreshold <- metadata(results)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)

        # Check that DEGs match the `DESeqResults` summary.
        match <- removeNA(str_match(
            string = capture.output(summary(results)),
            pattern = "^LFC.*\\s\\:\\s([0-9]+).*"
        ))
        assert_are_identical(
            x = length(up),
            y = as.integer(match[1L, 2L])
        )
        assert_are_identical(
            x = length(down),
            y = as.integer(match[2L, 2L])
        )

        TRUE
    }
)
