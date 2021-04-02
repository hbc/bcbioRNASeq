#' bcbio RNA-Seq data set
#'
#' `bcbioRNASeq` is an S4 class that extends `RangedSummarizedExperiment`, and
#' is designed to store a [bcbio](https://bcbio-nextgen.readthedocs.org) RNA-seq
#' analysis.
#'
#' @section Automatic metadata:
#'
#' The `metadata()` slot always contains:
#'
#' - Object version.
#' - bcbio data provenance information.
#' - File paths and timestamps.
#' - R session information.
#'
#' @author Michael Steinbaugh, Lorena Pantano
#' @note `bcbioRNASeq` extended `SummarizedExperiment` prior to v0.2.0, where we
#'   migrated to `RangedSummarizedExperiment`.
#' @note Updated 2021-04-02.
#' @export
setClass(
    Class = "bcbioRNASeq",
    contains = "RangedSummarizedExperiment"
)
setValidity(
    Class = "bcbioRNASeq",
    method = function(object) {
        metadata <- metadata(object)
        ## Changed from "version" to "packageVersion" on 2021-03-16.
        version <- metadata[["packageVersion"]]
        if (is.null(version)) {
            version <- metadata[["version"]]
        }
        ok <- validate(
            is(version, "package_version"),
            msg = "Package version is not defined in object."
        )
        if (!isTRUE(ok)) return(ok)
        ok <- validate(
            isTRUE(version >= 0.2),
            msg = "Object is older than v0.2, and cannot be easily updated."
        )
        if (!isTRUE(ok)) return(ok)
        ok <- validate(
            is(object, "RangedSummarizedExperiment"),
            hasDimnames(object),
            hasValidDimnames(object)
        )
        if (!isTRUE(ok)) return(ok)
        ## Metadata ------------------------------------------------------------
        ok <- validate(
            is.null(metadata[["metrics"]]),
            msg = "Legacy metrics detected inside object metadata."
        )
        ## Check that interesting groups defined in metadata are valid.
        ok <- validate(
            !isSubset(
                x = "interestingGroups",
                y = colnames(colData(object))
            ),
            msg = "'interestingGroups' column is not allowed in 'colData'."
        )
        if (!isTRUE(ok)) return(ok)
        ok <- validate(
            isSubset(
                x = metadata[["interestingGroups"]],
                y = colnames(colData(object))
            ),
            msg = "Interesting groups metadata not defined in 'colData'."
        )
        if (!isTRUE(ok)) return(ok)
        ## Error on legacy slot detection.
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
        ok <- validate(
            ## Using `as.logical()` here for R 3.4/BioC 3.6 compatibility.
            as.logical(!hasLength(intersect)),
            msg = sprintf("Legacy metadata: %s", toString(intersect))
        )
        if (!isTRUE(ok)) return(ok)
        ## Class checks (order independent).
        df <- "DataFrame"
        if (isTRUE(packageVersion("S4Vectors") >= "0.23")) {
            df <- c("DFrame", df)
        }
        ## Relaxed validity checks against `dataVersions`, `programsVersions`
        ## in v0.3.31, since this was reported to cause errors in upgrade from
        ## v0.2.9, which is still used by bcbio-nextgen conda recipe.
        ## See issue: https://github.com/hbc/bcbioRNASeq/issues/146
        ok <- validateClasses(
            object = metadata,
            expected = list(
                ## > "dataVersions" = df,
                ## Renamed "version" to "packageVersion" on 2021-03-16.
                ## > "packageVersion" = "package_version",
                ## > "programVersions" = df,
                "allSamples" = "logical",
                "bcbioCommandsLog" = "character",
                "bcbioLog" = "character",
                "call" = "call",
                "caller" = "character",
                "date" = "Date",
                "ensemblRelease" = "integer",
                "genomeBuild" = "character",
                "gffFile" = "character",
                "interestingGroups" = "character",
                "lanes" = "integer",
                "level" = "character",
                "organism" = "character",
                "projectDir" = "character",
                "runDate" = "Date",
                "sampleDirs" = "character",
                "sampleMetadataFile" = "character",
                "sessionInfo" = "session_info",
                "uploadDir" = "character",
                "wd" = "character",
                "yaml" = "list"
            ),
            subset = TRUE
        )
        if (!isTRUE(ok)) return(ok)
        ## Caller (and tximport) checks.
        ok <- validate(
            isSubset(metadata[["caller"]], .callers),
            isSubset(metadata[["level"]], .levels)
        )
        if (!isTRUE(ok)) return(ok)
        if (isSubset(metadata[["caller"]], .tximportCallers)) {
            ok <- validate(
                isSubset(
                    x = metadata[["countsFromAbundance"]],
                    y = eval(formals(tximport)[["countsFromAbundance"]])
                )
            )
            if (!isTRUE(ok)) return(ok)
            ok <- validateClasses(
                object = metadata,
                expected = list(
                    tx2gene = "Tx2Gene"
                ),
                subset = TRUE
            )
            if (!isTRUE(ok)) return(ok)
        }
        ## Assays --------------------------------------------------------------
        assayNames <- assayNames(object)
        ok <- validate(isSubset(.assays, assayNames))
        if (!isTRUE(ok)) return(ok)
        ## Check that all assays are matrices.
        ## Note that in previous versions, we slotted `DESeqDataSet` and
        ## `DESeqTransform`, which can result in metadata mismatches because
        ## those objects contain their own `colData` and `rowData`.
        ok <- validate(all(bapply(assays(object), is.matrix)))
        if (!isTRUE(ok)) return(ok)
        ## Caller-specific checks.
        caller <- metadata[["caller"]]
        ok <- validate(isString(caller))
        if (!isTRUE(ok)) return(ok)
        if (caller %in% .tximportCallers) {
            ok <- validate(isSubset(.tximportAssays, assayNames))
        } else if (caller %in% .featureCountsCallers) {
            ok <- validate(isSubset(.featureCountsAssays, assayNames))
        }
        if (!isTRUE(ok)) return(ok)
        ## Check for average transcript length matrix, if necessary.
        if (
            metadata[["caller"]] %in% .tximportCallers &&
            metadata[["countsFromAbundance"]] == "no"
        ) {
            ok <- validate(isSubset("avgTxLength", assayNames))
            if (!isTRUE(ok)) return(ok)
        }
        ## Row data ------------------------------------------------------------
        rowRanges <- rowRanges(object)
        rowData <- rowData(object)
        ok <- validate(
            is(rowRanges, "GRanges"),
            is(rowData, "DataFrame"),
            identical(names(rowRanges), rownames(object))
            ## This check fails on BioC 3.6; SummarizedExperiment 1.8.
            ## nolint start
            ## > identical(rownames(rowData), rownames(object))
            ## nolint end
        )
        if (!isTRUE(ok)) return(ok)
        ## Column data ---------------------------------------------------------
        colData <- colData(object)
        ok <- validate(
            identical(rownames(colData), colnames(object)),
            areDisjointSets(
                x = colnames(colData),
                y = c("sampleId", .legacyMetricsCols)
            )
        )
        if (!isTRUE(ok)) return(ok)
        ## Return.
        TRUE
    }
)
