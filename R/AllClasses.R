setClassUnion(
    name = "missingOrNULL",
    members = c("missing", "NULL")
)



# bcbioRNASeq ==================================================================
#' `bcbioRNASeq` Class
#'
#' `bcbioRNASeq` is an S4 class that extends `RangedSummarizedExperiment`, and
#' is designed to store a [bcbio](https://bcbio-nextgen.readthedocs.org) RNA-seq
#' analysis.
#'
#' @note `bcbioRNASeq` extended `SummarizedExperiment` prior to v0.2.0, where we
#'   migrated to `RangedSummarizedExperiment`.
#'
#' @family S4 Object
#' @author Michael Steinbaugh, Lorena Pantano
#' @export
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
#' x <- bcbioRNASeq(uploadDir)
#' print(x)
setClass(
    Class = "bcbioRNASeq",
    contains = "RangedSummarizedExperiment"
)



setValidity(
    Class = "bcbioRNASeq",
    method = function(object) {
        stopifnot(metadata(object)[["version"]] >= 0.2)
        assert_is_all_of(object, "RangedSummarizedExperiment")
        assert_has_dimnames(object)

        # Assays ---------------------------------------------------------------
        # Note that `rlog` and `vst` DESeqTransform objects are optional
        assert_is_subset(requiredAssays, assayNames(object))
        # Check that all assays are matrices
        assayCheck <- vapply(
            X = assays(object),
            FUN = is.matrix,
            FUN.VALUE = logical(1L),
            USE.NAMES = TRUE
        )
        if (!all(assayCheck)) {
            stop(paste(
                paste(
                    "Assays that are not matrix:",
                    toString(names(assayCheck[!assayCheck]))
                ),
                updateMessage,
                sep = "\n"
            ))
        }

        # Gene-level specific
        if (metadata(object)[["level"]] == "genes") {
            assert_is_subset("normalized", names(assays(object)))
        }

        # Row data -------------------------------------------------------------
        assert_is_all_of(rowRanges(object), "GRanges")
        assert_is_all_of(rowData(object), "DataFrame")

        # Column data ----------------------------------------------------------
        assert_are_disjoint_sets(colnames(colData(object)), legacyMetricsCols)

        # Metadata -------------------------------------------------------------
        metadata <- metadata(object)

        # Check that interesting groups defined in metadata are valid
        assert_is_subset(
            x = metadata[["interestingGroups"]],
            y = colnames(colData(object))
        )

        # Detect legacy metrics
        if (is.data.frame(metadata[["metrics"]])) {
            stop(paste(
                "`metrics` saved in `metadata()` instead of `colData()`.",
                updateMessage
            ))
        }

        # Detect legacy slots
        legacyMetadata <- c(
            "design",
            "ensemblVersion",
            "gtf",
            "gtfFile",
            "missingGenes",
            "programs",
            "yamlFile"
        )
        intersect <- intersect(names(metadata), legacyMetadata)
        if (length(intersect)) {
            stop(paste(
                paste(
                    "Legacy metadata slots:",
                    toString(sort(intersect))
                ),
                updateMessage,
                sep = "\n"
            ))
        }

        # v0.2.6: countsFromAbundance is now optional, since we're supporting
        # featureCounts aligned counts

        # Class checks (order independent)
        requiredMetadata <- list(
            allSamples = "logical",
            bcbioCommandsLog = "character",
            bcbioLog = "character",
            caller = "character",
            date = "Date",
            devtoolsSessionInfo = "session_info",
            ensemblRelease = "integer",
            genomeBuild = "character",
            gffFile = "character",
            interestingGroups = "character",
            lanes = "integer",
            level = "character",
            organism = "character",
            programVersions = "tbl_df",
            projectDir = "character",
            runDate = "Date",
            sampleDirs = "character",
            sampleMetadataFile = "character",
            template = "character",
            tx2gene = c("tx2gene", "data.frame"),
            uploadDir = "character",
            utilsSessionInfo = "sessionInfo",
            version = "package_version",
            wd = "character",
            yaml = "list"
        )
        classChecks <- invisible(mapply(
            name <- names(requiredMetadata),
            expected <- requiredMetadata,
            MoreArgs = list(metadata = metadata),
            FUN = function(name, expected, metadata) {
                actual <- class(metadata[[name]])
                if (!length(intersect(expected, actual))) {
                    FALSE
                } else {
                    TRUE
                }
            },
            SIMPLIFY = TRUE,
            USE.NAMES = TRUE
        ))
        if (!all(classChecks)) {
            stop(paste(
                "Metadata class checks failed.",
                updateMessage,
                printString(classChecks),
                sep = "\n"
            ))
        }

        # Additional assert checks
        # caller
        assert_is_subset(
            x = metadata[["caller"]],
            y = validCallers
        )
        # level
        assert_is_subset(
            x = metadata[["level"]],
            y = validLevels
        )
        # tx2gene
        tx2gene <- metadata[["tx2gene"]]
        # Switch to requiring `tx2gene` class in a future update.
        assertIsTx2gene(tx2gene)

        TRUE
    }
)



# DESeqAnalysis ================================================================
#' `DESeqAnalysis` Class
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
#' @family S4 Object
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @slot data `DESeqDataSet`.
#' @slot transform `DESeqTransform`.
#' @slot results `list`. One or more unshrunken `DESeqResults`.
#' @slot lfcShrink `list`. One or more shrunken `DESeqResults`.
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



#' `DESeqAnalysis` Generator
#'
#' @family S4 Object
#' @author Michael Steinbaugh
#' @export
#'
#' @param data `DESeqDataSet`.
#' @param transform `DESeqTransform`.
#' @param results `list`. One or more unshrunken `DESeqResults`. Assign the
#'   [DESeq2::results()] return here.
#' @param lfcShrink `list`. One or more shrunken `DESeqResults`. Assign the
#'   [DESeq2::lfcShrink()] return here.
#'
#' @return `DESeqAnalysis`.
#'
#' @examples
#' library(DESeq2)
#' dds <- as(bcb_small, "DESeqDataSet")
#' design(dds) <- ~ treatment
#' dds <- DESeq(dds)
#' class(dds)
#' vst <- varianceStabilizingTransformation(dds)
#' class(vst)
#' resultsNames(dds)
#' res <- results(dds, name = resultsNames(dds)[[2L]])
#' class(res)
#' x <- DESeqAnalysis(
#'     data = dds,
#'     transform = vst,
#'     results = list(res),
#'     lfcShrink = list(lfcShrink(dds = dds, coef = 2L))
#' )
#' print(x)
DESeqAnalysis <- function(
    data,
    transform,
    results,
    lfcShrink
) {
    new(
        Class = "DESeqAnalysis",
        data = data,
        transform = transform,
        results = results,
        lfcShrink = lfcShrink
    )
}



setValidity(
    Class = "DESeqAnalysis",
    method = function(object) {
        # Require valid dimnames for data set.
        assertHasValidDimnames(object@data)

        # Require gene-to-symbol mappings.
        assert_is_subset(
            x = c("geneID", "geneName"),
            y = colnames(rowData(object@data))
        )

        # DESeqDataSet and DESeqTransform must match.
        assert_are_identical(
            x = dimnames(object@data),
            y = dimnames(object@transform)
        )

        # DESeqDataSet and DESeqResults must match.
        invisible(lapply(
            X = object@results,
            FUN = function(res) {
                assert_are_identical(
                    x = rownames(res),
                    y = rownames(object@data)
                )
            }
        ))

        # DESeqResults and lfcShrink rownames must match, if shrinkage has been
        # calculated.
        if (length(object@lfcShrink) > 0L) {
            invisible(mapply(
                unshrunken = object@results,
                shrunken = object@lfcShrink,
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
#' `DESeqResultsTables` Class
#'
#' @family S4 Object
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#' @slot all `DESeqResults`. Original unmodified `DESeqResults`. Should contain
#'   all genes, including those with `NA` adjusted *P* values.
#' @slot deg `DataFrame`. Subset containing genes that pass adjusted *P* value
#'   and log2 fold change cutoffs.
#' @slot degUp `DataFrame`. Directional subset containing only upregulated
#'   genes.
#' @slot degDown `DataFrame`. Directional subset containing only downregulated
#'   genes.
#' @slot localFiles `list`. Local file paths.
#' @slot dropboxFiles `list`. Dropbox file paths.
setClass(
    Class = "DESeqResultsTables",
    slots = c(
        all = "DESeqResults",
        deg = "DataFrame",
        degUp = "DataFrame",
        degDown = "DataFrame",
        localFiles = "character",
        dropboxFiles = "list"
    ),
    # Consider setting an initialize method instead.
    prototype = list(
        localFiles = character(),
        dropboxFiles = list()
    )
)



#' `DESeqResultsTables` Generator
#'
#' @family S4 Object
#' @author Michael Steinbaugh
#' @export
#'
#' @param object `DESeqResults`.
#'
#' @return `DESeqResultsTables`.
#'
#' @examples
#' object <- deseq_small@results[[1L]]
#' class(object)
#' x <- DESeqResultsTables(object)
#' print(x)
DESeqResultsTables <- function(object) {
    resultsTables(object)
}



setValidity(
    Class = "DESeqResultsTables",
    method = function(object) {
        assert_is_a_string(contrastName(object@all))
        assert_is_a_number(metadata(object@all)[["alpha"]])
        assert_is_a_number(metadata(object@all)[["lfcThreshold"]])
        assert_is_subset(
            x = rownames(object@degUp),
            y = rownames(object@all)
        )
        assert_is_subset(
            x = rownames(object@degDown),
            y = rownames(object@all)
        )
        assert_are_disjoint_sets(
            x = rownames(object@degUp),
            y = rownames(object@degDown)
        )

        # Require that the subsets match the `DESeqResults` summary.
        match <- removeNA(str_match(
            string = capture.output(summary(object@all)),
            pattern = "^LFC.*\\s\\:\\s([0-9]+).*"
        ))
        assert_are_identical(
            x = nrow(object@degUp),
            y = as.integer(match[1, 2])
        )
        assert_are_identical(
            x = nrow(object@degDown),
            y = as.integer(match[2, 2])
        )

        # Object can have either local or Dropbox file paths, but not both.
        stopifnot(!all(
            length(object@localFiles) > 0L,
            length(object@dropboxFiles) > 0L
        ))

        TRUE
    }
)
