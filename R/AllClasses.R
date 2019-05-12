setClassUnion("missingOrNULL", c("missing", "NULL"))



#' @rdname bcbioRNASeq
#' @aliases NULL
#' @exportClass bcbioRNASeq
#' @usage NULL
bcbioRNASeq <- setClass(
    Class = "bcbioRNASeq",
    contains = "RangedSummarizedExperiment",
    validity = function(object) {
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
                bcbioBase::updateMessage,
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
                bcbioBase::updateMessage
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
                bcbioBase::updateMessage,
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
            tx2gene = "data.frame",
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
                bcbioBase::updateMessage,
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
        assertIsTx2gene(tx2gene)

        TRUE
    }
)
