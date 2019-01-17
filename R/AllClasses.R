setClassUnion(name = "missingOrNULL", members = c("missing", "NULL"))



#' bcbio RNA-Seq data set
#'
#' `bcbioRNASeq` is an S4 class that extends `RangedSummarizedExperiment`, and
#' is designed to store a [bcbio](https://bcbio-nextgen.readthedocs.org) RNA-seq
#' analysis.
#'
#' @section Automatic metadata:
#'
#' The [`metadata()`][S4Vectors::metadata] slot always contains:
#'
#' - Object version.
#' - bcbio data provenance information.
#' - File paths and timestamps.
#' - R session information.
#'
#' @author Michael Steinbaugh, Lorena Pantano, Rory Kirchner, Victor Barrera
#' @export
#'
#' @note `bcbioRNASeq` extended `SummarizedExperiment` prior to v0.2.0, where we
#'   migrated to `RangedSummarizedExperiment`.
setClass(
    Class = "bcbioRNASeq",
    contains = "RangedSummarizedExperiment",
    validity = function(object) {
        metadata <- metadata(object)

        # Return invalid for all objects older than v0.2.
        version <- metadata[["version"]]
        ok <- validate(
            is(version, "package_version"),
            version >= 0.2
        )
        if (!isTRUE(ok)) return(ok)

        ok <- validate(
            is(object, "RangedSummarizedExperiment"),
            hasDimnames(object)
        )
        if (!isTRUE(ok)) return(ok)

        # Metadata -------------------------------------------------------------
        ok <- validate(
            # Check for legacy metrics stashed in metadata, rather than defined
            # in colData.
            is.null(metadata[["metrics"]]),
            # Check that interesting groups defined in metadata are valid.
            isSubset(
                x = metadata[["interestingGroups"]],
                y = colnames(colData(object))
            )
        )
        if (!isTRUE(ok)) return(ok)

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
        ok <- validate(
            !hasLength(intersect),
            msg = sprintf("Legacy metadata: %s", toString(intersect))
        )
        if (!isTRUE(ok)) return(ok)

        # Class checks (order independent).
        ok <- validateClasses(
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
        if (!isTRUE(ok)) return(ok)

        # tximport checks.
        ok <- validate(
            isSubset(metadata[["caller"]], validCallers),
            isSubset(metadata[["level"]], validLevels),
            isSubset(
                x = metadata[["countsFromAbundance"]],
                y = eval(formals(tximport)[["countsFromAbundance"]])
            )
        )
        if (!isTRUE(ok)) return(ok)

        # Assays ---------------------------------------------------------------
        assayNames <- assayNames(object)
        ok <- validate(isSubset(requiredAssays, assayNames))
        if (!isTRUE(ok)) return(ok)

        # Check that all assays are matrices.
        # Note that in previous versions, we slotted `DESeqDataSet` and
        # `DESeqTransform`, which can result in metadata mismatches because
        # those objects contain their own `colData` and `rowData`.
        ok <- validate(all(bapply(assays(object), is.matrix)))
        if (!isTRUE(ok)) return(ok)

        # Caller-specific checks.
        caller <- metadata[["caller"]]
        ok <- validate(isString(caller))
        if (!isTRUE(ok)) return(ok)
        if (caller %in% tximportCallers) {
            ok <- validate(isSubset(tximportAssays, assayNames))
        } else if (caller %in% featureCountsCallers) {
            ok <- validate(isSubset(featureCountsAssays, assayNames))
        }
        if (!isTRUE(ok)) return(ok)

        # Check for average transcript length matrix, if necessary.
        if (
            metadata[["caller"]] %in% tximportCallers &&
            metadata[["countsFromAbundance"]] == "no"
        ) {
            ok <- validate(isSubset("avgTxLength", assayNames))
            if (!isTRUE(ok)) return(ok)
        }

        # Row data -------------------------------------------------------------
        rowRanges <- rowRanges(object)
        rowData <- rowData(object)
        ok <- validate(
            is(rowRanges, "GRanges"),
            is(rowData, "DataFrame")
        )
        if (!isTRUE(ok)) return(ok)

        if (hasLength(colnames(rowData))) {
            # Note that GTF/GFF annotations won't contain description.
            # The description column only gets returned via ensembldb.
            # This check will fail for bcbioRNASeq objects created prior to v0.3
            # update because we didn't use S4 Rle run-length encoding.
            ok <- validateClasses(
                object = rowData,
                expected = list(
                    broadClass = Rle,
                    geneBiotype = Rle,
                    geneID = Rle,
                    geneName = Rle
                ),
                subset = TRUE
            )
            if (!isTRUE(ok)) return(ok)
        }

        # Column data ----------------------------------------------------------
        colData <- colData(object)
        ok <- validate(
            isSubset("sampleName", colnames(colData)),
            areDisjointSets(colnames(colData), legacyMetricsCols)
        )
        if (!isTRUE(ok)) return(ok)

        TRUE
    }
)
