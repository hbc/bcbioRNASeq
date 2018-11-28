setClassUnion(name = "missingOrNULL", members = c("missing", "NULL"))



.valid <- function(list) {
    invalid <- Filter(f = Negate(isTRUE), x = list)
    if (has_length(invalid)) {
        unlist(invalid)
    } else {
        TRUE
    }
}



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
#' The `metadata()` slot always contains:
#'
#' - Object version.
#' - bcbio data provenance information.
#' - File paths and timestamps.
#' - R session information.
#'
#' @family S4 classes
#' @author Michael Steinbaugh, Lorena Pantano
#' @export
#'
#' @seealso `bcbioRNASeq()`.
setClass(Class = "bcbioRNASeq", contains = "RangedSummarizedExperiment")
setValidity(
    Class = "bcbioRNASeq",
    method = function(object) {
        valid <- list()
        metadata <- metadata(object)

        # Return invalid for all objects older than v0.2.
        version <- metadata[["version"]]
        valid[["version"]] <- validate_that(
            is(version, "package_version"),
            version >= 0.2
        )

        valid[["rse"]] <- validate_that(
            is(object, "RangedSummarizedExperiment")
        )

        valid[["dimnames"]] <- validate_that(
            has_dimnames(object)
        )

        # Metadata -------------------------------------------------------------
        # Check for legacy metrics.
        valid[["metrics"]] <- validate_that(
            !is.data.frame(metadata[["metrics"]]),
            msg = "`metrics` saved in `metadata()` instead of `colData()`."
        )

        # Check that interesting groups defined in metadata are valid.
        valid[["interestingGroups"]] <- validate_that(
            is_subset(
                x = metadata[["interestingGroups"]],
                y = colnames(colData(object))
            )
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
        valid[["legacyMetadata"]] <- validate_that(
            !has_length(intersect),
            msg = paste(
                "Legacy metadata slots:",
                toString(sort(intersect)),
                sep = "\n"
            )
        )

        # Class checks (order independent).
        valid[["metadata"]] <- validateClasses(
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
        valid[["metadata2"]] <- validate_that(
            is_subset(metadata[["caller"]], validCallers),
            is_subset(metadata[["level"]], validLevels)
        )

        if (is.character(metadata[["countsFromAbundance"]])) {
            valid[["countsFromAbundance"]] <- validate_that(
                is_subset(
                    x = metadata[["countsFromAbundance"]],
                    y = eval(formals(tximport)[["countsFromAbundance"]])
                )
            )
        }

        # Assays ---------------------------------------------------------------
        assayNames <- assayNames(object)
        valid[["assayNames"]] <- validate_that(
            is_subset(requiredAssays, assayNames)
        )

        # Check that all assays are matrices.
        # Note that in previous versions, we slotted `DESeqDataSet` and
        # `DESeqTransform`, which can result in metadata mismatches because
        # those objects contain their own `colData()` and `rowData()`.
        isMatrix <- vapply(
            X = assays(object),
            FUN = is.matrix,
            FUN.VALUE = logical(1L),
            USE.NAMES = TRUE
        )
        valid[["assays"]] <- validate_that(
            all(isMatrix),
            msg = paste(
                "Assays that are not matrix:",
                toString(names(valid[!valid])),
                sep = "\n"
            )
        )

        # Caller-specific checks.
        caller <- metadata[["caller"]]
        if (caller %in% tximportCallers) {
            valid[["caller"]] <- validate_that(
                is_subset(tximportAssays, assayNames)
            )
        } else if (caller %in% featureCountsCallers) {
            valid[["caller"]] <- validate_that(
                is_subset(featureCountsAssays, assayNames)
            )
        }

        # Check for average transcript length matrix, if necessary.
        if (
            metadata[["caller"]] %in% tximportCallers &&
            metadata[["countsFromAbundance"]] == "no"
        ) {
            valid[["avgTxLength"]] <- validate_that(
                is_subset("avgTxLength", assayNames)
            )
        }

        # Row data -------------------------------------------------------------
        valid[["rowRanges"]] <- validate_that(
            is(rowRanges(object), "GRanges")
        )
        rowData <- rowData(object)
        if (has_length(colnames(rowData))) {
            # Note that GTF/GFF annotations won't contain description.
            valid[["rowData"]] <- validateClasses(
                object = rowData,
                expected = list(
                    broadClass = Rle,
                    geneBiotype = Rle,
                    geneID = Rle,
                    geneName = Rle
                ),
                subset = TRUE
            )
        }

        # Column data ----------------------------------------------------------
        colData <- colData(object)
        valid[["colData"]] <- validate_that(
            is_subset("sampleName", colnames(colData)),
            are_disjoint_sets(colnames(colData), legacyMetricsCols)
        )

        .valid(list = valid)
    }
)
