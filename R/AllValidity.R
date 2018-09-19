# bcbioRNASeq ==================================================================
.setValidity.bcbioRNASeq <-  # nolint
    function(object) {
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

        # Support for legacy `devtoolsSessionInfo` stash, which has been
        # renamed to simply `sessionInfo`. Previously, we stashed both
        # `devtoolsSessionInfo` and `utilsSessionInfo`.
        if ("devtoolsSessionInfo" %in% names(metadata)) {
            metadata[["sessionInfo"]] <- metadata[["devtoolsSessionInfo"]]
        }

        # Check that interesting groups defined in metadata are valid.
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

        # Detect legacy slots.
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
        # featureCounts aligned counts.

        # Class checks (order independent)
        requiredMetadata <- list(
            allSamples = "logical",
            bcbioCommandsLog = "character",
            bcbioLog = "character",
            caller = "character",
            date = "Date",
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
            sessionInfo = "session_info",
            uploadDir = "character",
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



setValidity(
    Class = "bcbioRNASeq",
    method = .setValidity.bcbioRNASeq
)



# DESeqAnalysis ================================================================
.setValidity.DESeqAnalysis <-  # nolint
    function(object) {
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



setValidity(
    Class = "DESeqAnalysis",
    method = .setValidity.DESeqAnalysis
)



# DESeqResultsTables ===========================================================
.setValidity.DESeqResultsTables <-  # nolint
    function(object) {
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



setValidity(
    Class = "DESeqResultsTables",
    method = .setValidity.DESeqResultsTables
)
