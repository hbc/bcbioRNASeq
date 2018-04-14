#' Update Object
#'
#' Update old objects created by the bcbioRNASeq package. The session
#' information metadata is preserved from the time when the bcbio data was
#' originally loaded into R.
#'
#' @section Legacy `bcbioRNADataSet` class:
#' Support for `bcbioRNADataSet` objects was dropped in v0.2.0 of the package.
#' If you need to load one of these objects, please install an older release.
#'
#' @section Legacy objects created with `bcbioRnaseq`:
#' The previous `bcbioRnaseq` package (note case) must be reinstalled to load
#' objects from versions <= 0.0.20. We changed the name of the package to
#' `bcbioRNASeq` starting in v0.0.21.
#'
#' @name updateObject
#' @family S4 Class Definition
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics updateObject
#'
#' @inheritParams general
#' @param rowRanges `GRanges` object that defines the row annotations. Since
#'   we converted to `RangedSummarizedExperiment` in v0.2.0, this option had
#'   to be added to enable updating of newly required [rowRanges()] slot.
#'
#' @return `bcbioRNASeq`.
#'
#' @examples
#' updateObject(bcb_small)
NULL



# Methods ======================================================================
#' @rdname updateObject
#' @export
setMethod(
    "updateObject",
    signature("bcbioRNASeq"),
    function(
        object,
        rowRanges
    ) {
        version <- metadata(object)[["version"]]
        assert_is_all_of(version, c("package_version", "numeric_version"))
        message(paste("Upgrading from", version, "to", packageVersion))

        # Check for legacy bcbio slot
        if (.hasSlot(object, "bcbio")) {
            message("Legacy bcbio slot detected")
        }

        # Coerce to RangedSummarizedExperiment
        rse <- as(object, "RangedSummarizedExperiment")

        # Assays ===============================================================
        # Using `slot()` here to avoid error on missing rowRanges
        assays <- slot(rse, "assays")
        if (!all(assayNames(rse) %in% requiredAssays)) {
            # Coerce assays to a standard list
            assays <- lapply(seq_along(assays), function(a) {
                assays[[a]]
            })
            names(assays) <- assayNames(rse)
            # DESeqDataSet (only on develop branch, safe to remove later)
            if ("dds" %in% names(assays)) {
                message("Dropping legacy DESeqDataSet from assays")
                assays[["dds"]] <- NULL
            }
            # length (from tximport)
            if (!"length" %in% names(assays)) {
                message("Moving length matrix from legacy bcbio slot to assays")
                length <- slot(object, "bcbio")[["tximport"]][["length"]]
                assert_is_matrix(length)
                assays[["length"]] <- length
            }
            # DESeq2 normalized counts
            if ("normalized" %in% names(assays)) {
                message("Dropping legacy DESeq2 normalized counts from assays")
                assays[["normalized"]] <- NULL
            }
            # rlog
            if (is(assays[["rlog"]], "DESeqTransform")) {
                message("Coercing rlog DESeqTransform to matrix in assays")
                assays[["rlog"]] <- assay(assays[["rlog"]])
            }
            # tmm
            if ("tmm" %in% names(assays)) {
                message("Dropping tmm from assays. Now calculated on the fly.")
                assays[["tmm"]] <- NULL
            }
            # vst
            if (is(assays[["vst"]], "DESeqTransform")) {
                message("Coercing vst DESeqTransform to matrix in assays")
                assays[["vst"]] <- assay(assays[["vst"]])
            }
            assays <- Filter(Negate(is.null), assays)
            # Put the required assays first, in order
            assays <- assays[unique(c(requiredAssays, names(assays)))]
            assert_is_subset(requiredAssays, names(assays))
        }

        # Row data =============================================================
        if (.hasSlot(rse, "rowRanges")) {
            rowRanges <- slot(rse, "rowRanges")
        }
        assert_is_all_of(rowRanges, "GRanges")

        # Column data ==========================================================
        colData <- colData(rse)
        colData <- sanitizeSampleData(colData)

        # Metadata =============================================================
        metadata <- metadata(rse)

        # bcbioLog
        if (is.null(metadata[["bcbioLog"]])) {
            message("Setting bcbioLog as empty character")
            metadata[["bcbioLog"]] <- character()
        }

        # bcbioCommandsLog
        if (is.null(metadata[["bcbioCommandsLog"]])) {
            message("Setting bcbioCommands as empty character")
            metadata[["bcbioCommandsLog"]] <- character()
        }

        # caller
        if (!"caller" %in% names(metadata)) {
            message("Setting caller as salmon")
            metadata[["caller"]] <- "salmon"
        }

        # countsFromAbundance
        if (!"countsFromAbundance" %in% names(metadata)) {
            message("Setting countsFromAbundance as lengthScaledTPM")
            metadata[["countsFromAbundance"]] <- "lengthScaledTPM"
        }

        # design
        if ("design" %in% names(metadata)) {
            message("Dropping legacy design formula")
            metadata[["design"]] <- NULL
        }

        # ensemblRelease
        if ("ensemblVersion" %in% names(metadata)) {
            # Renamed in v0.2.0
            message("Renaming ensemblVersion to ensemblRelease")
            metadata[["ensemblRelease"]] <- metadata[["ensemblVersion"]]
            metadata[["ensemblVersion"]] <- NULL
        }
        if (!is.integer(metadata[["ensemblRelease"]])) {
            message("Setting ensemblRelease as integer")
            metadata[["ensemblRelease"]] <-
                as.integer(metadata[["ensemblRelease"]])
        }

        # genomeBuild
        if (!is.character(metadata[["genomeBuild"]])) {
            message("Setting genomeBuild as empty character")
            metadata[["genomeBuild"]] <- character()
        }

        # gffFile
        if ("gtfFile" %in% names(metadata)) {
            message("Renaming gtfFile to gffFile")
            metadata[["gffFile"]] <- metadata[["gtfFile"]]
            metadata[["gtfFile"]] <- NULL
        }
        if (!"gffFile" %in% names(metadata)) {
            message("Setting gffFile as empty character")
            metadata[["gffFile"]] <- character()
        }

        # gtf
        if ("gtf" %in% names(metadata)) {
            message("Removing stashed GTF")
            metadata <- metadata[setdiff(names(metadata), "gtf")]
        }

        # lanes
        if (!is.integer(metadata[["lanes"]])) {
            message("Setting lanes as integer")
            metadata[["lanes"]] <- as.integer(metadata[["lanes"]])
        }

        # level
        if (!"level" %in% names(metadata)) {
            message("Setting level as genes")
            metadata[["level"]] <- "genes"
        }

        # metrics
        if (length(intersect(
            colnames(metadata[["metrics"]]), legacyMetricsCols
        ))) {
            metrics <- metadata[["metrics"]]
            assert_is_data.frame(metrics)

            # Drop sample name columns
            legacyNameCols <- c(bcbioBase::metadataPriorityCols, "name")
            if (length(intersect(colnames(metrics), legacyNameCols))) {
                message("Dropping legacy sample names from metrics")
                metrics <- metrics %>%
                    .[, setdiff(colnames(.), legacyNameCols), drop = FALSE]
            }

            # Rename 5'3' bias
            if ("x53Bias" %in% colnames(metrics)) {
                message("Renaming x53Bias to x5x3Bias")
                metrics[["x5x3Bias"]] <- metrics[["x53Bias"]]
                metrics[["x53Bias"]] <- NULL
            }

            # Rename rRNA rate
            if (!"rrnaRate" %in% colnames(metrics)) {
                col <- grep(
                    pattern = "rrnarate",
                    x = colnames(metrics),
                    ignore.case = TRUE,
                    value = TRUE
                )
                assert_is_a_string(col)
                message(paste("Renaming", col, "to rrnaRate"))
                metrics[["rrnaRate"]] <- metrics[[col]]
                metrics[[col]] <- NULL
            }

            metadata[["metrics"]] <- metrics
        }

        # programVersions
        if (!"programVersions" %in% names(metadata) &&
            "programs" %in% names(metadata)) {
            message("Renaming programs to programVersions")
            metadata[["programVersions"]] <- metadata[["programs"]]
            metadata <- metadata[setdiff(names(metadata), "programs")]
        }

        # rowRangesMetadata
        if (!"rowRangesMetadata" %in% names(metadata)) {
            message("Setting rowRangesMetadata as empty tibble")
            metadata[["rowRangesMetadata"]] <- tibble()
        }

        # sampleMetadataFile
        if (!is.character(metadata[["sampleMetadataFile"]])) {
            message("Setting sampleMetadataFile as empty character")
            metadata[["sampleMetadataFile"]] <- character()
        }

        # tx2gene
        if (any(c("enstxp", "ensgene") %in% colnames(metadata[["tx2gene"]]))) {
            message("Renaming enstxp, ensgene to txID, geneID")
            x <- metadata[["tx2gene"]]
            x[["txID"]] <- x[["enstxp"]]
            x[["enstxp"]] <- NULL
            x[["geneID"]] <- x[["ensgene"]]
            x[["ensgene"]] <- NULL
            metadata[["tx2gene"]] <- x
        }

        # unannotatedGenes
        if ("missingGenes" %in% names(metadata)) {
            message("Renaming missingGenes to unannotatedGenes")
            metadata[["unannotatedGenes"]] <- metadata[["missingGenes"]]
            metadata[["missingGenes"]] <- NULL
        }

        # yamlFile
        if ("yamlFile" %in% names(metadata)) {
            message("Dropping yamlFile path")
            metadata[["yamlFile"]] <- NULL
        }

        # version
        metadata[["previousVersion"]] <- metadata[["version"]]
        metadata[["version"]] <- packageVersion

        # Return ===============================================================
        .new.bcbioRNASeq(
            assays = assays,
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
    }
)
