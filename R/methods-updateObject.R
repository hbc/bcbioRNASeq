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
#' loadRemoteData("http://bcbiornaseq.seq.cloud/f1000v1/bcb.rda")
#' metadata(bcb)[["version"]]
#' organism <- metadata(bcb)[["organism"]]
#' rowRanges <- genes(organism)
#' updateObject(bcb, rowRanges = rowRanges)
NULL



# Constructors =================================================================
.updateObject.bcbioRNASeq <- function(  # nolint
    object,
    rowRanges
) {
    version <- metadata(object)[["version"]]
    assert_is_all_of(version, c("package_version", "numeric_version"))
    inform(paste("Upgrading from", version, "to", packageVersion))

    # RangedSummarizedExperiment
    rse <- as(object, "RangedSummarizedExperiment")
    assays <- slot(rse, "assays")
    if (.hasSlot(rse, "rowRanges")) {
        rowRanges <- slot(rse, "rowRanges")
    }
    assert_is_all_of(rowRanges, "GRanges")
    colData <- colData(rse)
    metadata <- metadata(rse)

    # Assays ===================================================================
    if (!all(assayNames(rse) %in% requiredAssays)) {
        # Ensure assays is a standard list
        assays <- lapply(seq_along(assays), function(a) {
            assays[[a]]
        })
        names(assays) <- assayNames(rse)
        if (!"dds" %in% assayNames(rse)) {
            # DESeqDataSet
            dds <- slot(object, "bcbio")[["DESeqDataSet"]]
            assert_is_all_of(dds, "DESeqDataSet")
            validObject(dds)
            assays[["dds"]] <- dds
        }
        if (!"length" %in% assayNames(rse)) {
            # tximport length
            length <- slot(object, "bcbio")[["tximport"]][["length"]]
            assert_is_matrix(length)
            assays[["length"]] <- length
        }
        if ("normalized" %in% assayNames(rse)) {
            # Legacy DESeq2 normalized counts
            inform("Dropping normalized from assays")
            assays[["normalized"]] <- NULL
        }
        if ("tmm" %in% assayNames(rse)) {
            # Calculate TMM on the fly instead
            inform("Dropping tmm from assays")
            assays[["tmm"]] <- NULL
        }
        assays <- Filter(Negate(is.null), assays)
        assays <- assays[unique(c(requiredAssays, names(assays)))]
        assert_is_subset(requiredAssays, names(assays))
    }

    # Metadata upgrades ========================================================
    metadata <- metadata(rse)

    # bcbioLog
    if (is.null(metadata[["bcbioLog"]])) {
        inform("Setting bcbioLog as empty character")
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
    if (!"design" %in% names(metadata)) {
        # This should only be stashed inside DESeqDataSet.
        # Formulas inside the metadata slot blow up memory and disk usage in R.
        inform("Dropping design formula inside metadata")
        metadata[["design"]] <- NULL
    }

    # ensemblRelease
    if ("ensemblVersion" %in% names(metadata)) {
        # Renamed in v0.2.0
        inform("Renaming ensemblVersion to ensemblRelease")
        metadata[["ensemblRelease"]] <- metadata[["ensemblVersion"]]
        metadata[["ensemblVersion"]] <- NULL
    }
    if (!is.integer(metadata[["ensemblRelease"]])) {
        inform("Setting ensemblRelease as integer")
        metadata[["ensemblRelease"]] <- as.integer(metadata[["ensemblRelease"]])
    }

    # genomeBuild
    if (!is.character(metadata[["genomeBuild"]])) {
        inform("Setting genomeBuild as empty character")
        metadata[["genomeBuild"]] <- character()
    }

    # gffFile
    if (!"gffFile" %in% names(metadata)) {
        inform("Setting gffFile as empty character")
        metadata[["gffFile"]] <- character()
    }
    if (!"gtfFile" %in% names(metadata)) {
        inform("Setting gffFile as gtfFile")
        metadata[["gffFile"]] <- metadata[["gtfFile"]]
        metadata[["gtfFile"]] <- NULL
    }

    # gtf
    if ("gtf" %in% names(metadata)) {
        inform("Removing stashed GTF")
        metadata <- metadata[setdiff(names(metadata), "gtf")]
    }

    # isSpike
    if (!is.character(metadata[["isSpike"]])) {
        inform("Setting isSpike as empty character")
        metadata[["isSpike"]] <- character()
    }

    # lanes
    if (!is.integer(metadata[["lanes"]])) {
        inform("Setting lanes as integer")
        metadata[["lanes"]] <- as.integer(metadata[["lanes"]])
    }

    # level
    if (!"level" %in% names(metadata)) {
        inform("Setting level as genes")
        metadata[["level"]] <- "genes"
    }

    # metrics
    if (length(intersect(colnames(metadata[["metrics"]]), legacyMetricsCols))) {
        metrics <- metadata[["metrics"]]
        assert_is_data.frame(metrics)

        # Drop sample name columns
        legacyNameCols <- c(metadataPriorityCols, "name")
        if (length(intersect(colnames(metrics), legacyNameCols))) {
            inform("Dropping legacy sample names from metrics")
            metrics <- metrics %>%
                .[, setdiff(colnames(.), legacyNameCols), drop = FALSE]
        }

        # Rename 5'3' bias
        if ("x53Bias" %in% colnames(metrics)) {
            inform("Renaming x53Bias to x5x3Bias")
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
            inform(paste("Renaming", col, "to rrnaRate"))
            metrics[["rrnaRate"]] <- metrics[[col]]
            metrics[[col]] <- NULL
        }

        metadata[["metrics"]] <- metrics
    }

    # programVersions
    if (!"programVersions" %in% names(metadata) &&
        "programs" %in% names(metadata)) {
        inform("Renaming programs to programVersions")
        metadata[["programVersions"]] <- metadata[["programs"]]
        metadata <- metadata[setdiff(names(metadata), "programs")]
    }

    # rowRangesMetadata
    if (!"rowRangesMetadata" %in% names(metadata)) {
        inform("Setting rowRangesMetadata as empty tibble")
        metadata[["rowRangesMetadata"]] <- tibble()
    }

    # sampleMetadataFile
    if (!is.character(metadata[["sampleMetadataFile"]])) {
        inform("Setting sampleMetadataFile as empty character")
        metadata[["sampleMetadataFile"]] <- character()
    }

    # transformationLimit
    if (!is.numeric(metadata[["transformationLimit"]])) {
        inform("Setting transformationLimit as Inf")
        metadata[["transformationLimit"]] <- Inf
    }

    # unannotatedGenes
    if ("missingGenes" %in% names(metadata)) {
        inform("Renaming missingGenes to unannotatedGenes")
        metadata[["unannotatedGenes"]] <- metadata[["missingGenes"]]
        metadata[["missingGenes"]] <- NULL
    }

    # yamlFile
    if ("yamlFile" %in% names(metadata)) {
        inform("Dropping yamlFile path")
        metadata[["yamlFile"]] <- NULL
    }

    metadata[["version"]] <- packageVersion
    metadata[["previousVersion"]] <- metadata[["version"]]
    metadata[["upgradeDate"]] <- Sys.Date()

    # Drop any NULL legacy metadata slots
    metadata <- Filter(Negate(is.null), metadata)

    # Return ===================================================================
    # Generate RangedSummarizedExperiment. Resizes rowRanges dynamically.
    rse <- prepareSummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata
    )
    validObject(rse)

    # Return ===================================================================
    to <- new("bcbioRNASeq", rse)
    validObject(to)
    to
}



# Methods ======================================================================
#' @rdname updateObject
#' @export
setMethod(
    "updateObject",
    signature("bcbioRNASeq"),
    .updateObject.bcbioRNASeq
)
