# TODO Consider a utility that we can use in `[` and `updateObject`.
# .updateMetadata <- function(metadata) {
#     metadata[["template"]] <- NULL
# }



#' Update an Object to Its Current Class Definition
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
#' @family S4 Functions
#' @author Michael Steinbaugh
#' @importFrom BiocGenerics updateObject
#' @export
#'
#' @inheritParams general
#' @param rowRanges `GRanges` or `NULL`. Row annotations. Since we converted to
#'   `RangedSummarizedExperiment` in v0.2.0, this option had to be added to
#'   enable updating of newly required [rowRanges()] slot. Objects that are >=
#'   v0.2 don't require this argument and it can be left `NULL`.
#'
#' @return `bcbioRNASeq`.
#'
#' @examples
#' data(bcb_small)
#' updateObject(bcb_small)
NULL



.updateObject.bcbioRNASeq <-  # nolint
    function(
        object,
        rowRanges = NULL
    ) {
        version <- slot(object, "metadata")[["version"]]
        assert_is_all_of(version, c("package_version", "numeric_version"))
        message(paste0(
            "Upgrading from ", version, " to ", packageVersion, "..."
        ))

        # Check for legacy bcbio slot.
        if (.hasSlot(object, "bcbio")) {
            message("Legacy bcbio slot detected.")
        }

        # Assays ---------------------------------------------------------------
        # Coerce the assays (e.g. ShallowSimpleListAssays) back to list.
        # Using `slot()` here to avoid error on missing rowRanges.
        # Ensure this comes before the rowRanges handling (see below).
        assays <- slot(object, "assays")
        assays <- lapply(seq_along(assays), function(a) {
            assays[[a]]
        })
        names(assays) <- assayNames(object)

        caller <- metadata(object)[["caller"]]
        assert_is_a_string(caller)
        assert_is_subset(caller, validCallers)

        level <- metadata(object)[["level"]]
        assert_is_a_string(level)
        assert_is_subset(level, validLevels)

        # Rename main assay from "raw" to "counts".
        if ("raw" %in% names(assays)) {
            message("Renaming primary assay to counts.")
            assays[["counts"]] <- assays[["raw"]]
            assays[["raw"]] <- NULL
        }

        # Legacy average transcript length matrix.
        countsFromAbundance <- metadata(object)[["countsFromAbundance"]]
        if (
            "length" %in% names(assays) &&
            countsFromAbundance != "no"
        ) {
            # Note that if we use `countsFromAbundance = "no"`, this will be
            # imported as `avgTxLength`, which integrates with DESeq2. We don't
            # want to keep average transcript lengths in assays when using
            # `countsFromAbundance = "lengthScaledTPM"`, since these have
            # already been applied to the raw counts matrix.
            message(paste0(
                "Dropping average transcript length matrix from assays.\n",
                "Counts were imported with tximport using ",
                "`countsFromAbundance = ", deparse(countsFromAbundance), "`. ",
                "Since the counts are already scaled, ",
                "inclusion of the length matrix isn't recommended."
            ))
            assays[["length"]] <- NULL
        }
        rm(countsFromAbundance)

        # Legacy TMM counts.
        if ("tmm" %in% names(assays)) {
            message(paste(
                "Dropping tmm from assays.",
                "Now calculated on the fly instead."
            ))
            assays[["tmm"]] <- NULL
        }

        # Gene-level-specific assays (DESeq2).
        if (level == "genes") {
            # DESeq2 normalized counts.
            if (is(assays[["normalized"]], "DESeqDataSet")) {
                dds <- assays[["normalized"]]
                assays[["normalized"]] <- counts(dds, normalized = TRUE)
            } else if (!"normalized" %in% names(assays)) {
                dds <- .new.DESeqDataSetFromMatrix(assays[["counts"]])
                dds <- DESeq(dds)
                assays[["normalized"]] <- counts(dds, normalized = TRUE)
            }

            # Variance-stabilizing transformation.
            if (is(assays[["vst"]], "DESeqTransform")) {
                message("Coercing vst DESeqTransform to matrix in assays.")
                assays[["vst"]] <- assay(assays[["vst"]])
            }

            # Regularized log.
            if (is(assays[["rlog"]], "DESeqTransform")) {
                message("Coercing rlog DESeqTransform to matrix in assays.")
                assays[["rlog"]] <- assay(assays[["rlog"]])
            }
        }

        # Always put the required assays first.
        assays <- assays[unique(c(requiredAssays, names(assays)))]
        assert_is_subset(requiredAssays, names(assays))

        # Check for required caller-specific assays.
        if (caller %in% tximportCallers) {
            assert_is_subset(tximportAssays, names(assays))
        } else if (caller %in% featureCountsCallers) {
            assert_is_subset(featureCountsAssays, names(assays))
        }

        # Row data -------------------------------------------------------------
        # Error if the object already has rowRanges, and the user is trying
        # to supply new ones.
        if (
            .hasSlot(object, "rowRanges") &&
            !is.null(rowRanges)
        ) {
            stop(paste(
                "Object already contains rowRanges.",
                "Don't attempt to slot new ones with `rowRanges` argument.",
                sep = "\n"
            ), call. = FALSE)
        } else if (.hasSlot(object, "rowRanges")) {
            rowRanges <- rowRanges(object)
        } else if (is.null(rowRanges)) {
            warning(paste(
                "`rowRanges` are now recommended for gene annotations."
            ), call. = FALSE)
            # Generate empty ranges if not supplied by the user.
            rowRanges <- emptyRanges(names = rownames(assays[[1L]]))
        }
        assert_is_all_of(rowRanges, "GRanges")

        # Column data ----------------------------------------------------------
        colData <- colData(object)

        # Move metrics from metadata into colData, if necessary.
        metrics <- metadata(object)[["metrics"]]
        if (!is.null(metrics)) {
            message("Moving metrics from `metadata()` into `colData()`.")
            assert_is_data.frame(metrics)

            # Always remove legacy `name` column.
            metrics[["name"]] <- NULL

            # Rename 5'3' bias.
            if ("x53Bias" %in% colnames(metrics)) {
                message("Renaming x53Bias to x5x3Bias.")
                metrics[["x5x3Bias"]] <- metrics[["x53Bias"]]
                metrics[["x53Bias"]] <- NULL
            }

            # Rename rRNA rate.
            if (!"rrnaRate" %in% colnames(metrics)) {
                col <- grep(
                    pattern = "rrnarate",
                    x = colnames(metrics),
                    ignore.case = TRUE,
                    value = TRUE
                )
                assert_is_a_string(col)
                message(paste("Renaming", col, "to rrnaRate."))
                metrics[["rrnaRate"]] <- metrics[[col]]
                metrics[[col]] <- NULL
            }

            # Only include columns not already present in colData.
            setdiff <- setdiff(colnames(metrics), colnames(colData))
            metrics <- metrics[, sort(setdiff), drop = FALSE]

            colData <- cbind(colData, metrics)
            metadata(object)[["metrics"]] <- NULL
        }

        # Remove legacy `sampleID` and `description` columns, if present.
        colData[["sampleID"]] <- NULL
        colData[["description"]] <- NULL

        # Metadata -------------------------------------------------------------
        metadata <- metadata(object)

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
        if (!is(metadata[["tx2gene"]], "Tx2Gene")) {
            metadata[["tx2gene"]] <- tx2gene(metadata[["tx2gene"]])
        }

        # Dead genes: "missing" or "unannotated"
        if ("missingGenes" %in% names(metadata)) {
            message("Dropping missingGenes from metadata.")
            metadata[["missingGenes"]] <- NULL
        }
        if ("unannotatedGenes" %in% names(metadata)) {
            message("Dropping unannotatedGenes from metadata.")
            metadata[["unannotatedGenes"]] <- NULL
        }

        # yamlFile
        if ("yamlFile" %in% names(metadata)) {
            message("Dropping yamlFile path.")
            metadata[["yamlFile"]] <- NULL
        }

        # version
        metadata[["previousVersion"]] <- metadata[["version"]]
        metadata[["version"]] <- packageVersion

        # Return ---------------------------------------------------------------
        .new.bcbioRNASeq(
            assays = assays,
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
    }



#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature("bcbioRNASeq"),
    definition = .updateObject.bcbioRNASeq
)
