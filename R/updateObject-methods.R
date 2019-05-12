#' Update an object to its current class definition
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
#' @family S4 Object
#' @author Michael Steinbaugh
#'
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
#' updateObject(bcb_small)
NULL



updateObject.bcbioRNASeq <-  # nolint
    function(
        object,
        rowRanges = NULL
    ) {
        version <- slot(object, "metadata")[["version"]]
        assert_is_all_of(version, c("package_version", "numeric_version"))
        message(paste("Upgrading from", version, "to", packageVersion))

        # Check for legacy bcbio slot
        if (.hasSlot(object, "bcbio")) {
            message("Legacy bcbio slot detected")
        }

        # Coerce assays to list ------------------------------------------------
        # Ensure this comes before the rowRanges handling
        # Using `slot()` here to avoid error on missing rowRanges
        assays <- slot(object, "assays")
        # Coerce ShallowSimpleListAssays to list
        assays <- lapply(seq_along(assays), function(a) {
            assays[[a]]
        })
        names(assays) <- assayNames(object)

        # Row data -------------------------------------------------------------
        rownames <- rownames(assays[[1L]])
        # This section needs to come before the assay modifications
        if (.hasSlot(object, "rowRanges")) {
            hasRowRanges <- TRUE
            intRowRanges <- slot(object, "rowRanges")
        } else {
            hasRowRanges <- FALSE
            # Generate empty genomic ranges if not supplied by the user
            intRowRanges <- emptyRanges(names = rownames)
        }
        assert_is_all_of(intRowRanges, "GRanges")

        # Regenerate RangedSummarizedExperiment --------------------------------
        rse <- SummarizedExperiment(
            assays = assays,
            rowRanges = intRowRanges,
            colData = colData(object),
            metadata = metadata(object)
        )
        validObject(rse)
        rm(assays, intRowRanges)
        if (is.null(rowRanges)) {
            if (!isTRUE(hasRowRanges)) {
                warning("`rowRanges` are now recommended for gene annotations")
            }
            rowRanges <- rowRanges(rse)
        }

        # Metadata -------------------------------------------------------------
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
        if ("txID" %in% colnames(metadata[["tx2gene"]])) {
            message("tx2gene: Renaming `txID` to `transcriptID`")
            assert_are_identical(
                c("txID", "geneID"),
                colnames(metadata[["tx2gene"]])
            )
            colnames(metadata[["tx2gene"]]) <- c("transcriptID", "geneID")
        }
        if (any(c("enstxp", "ensgene") %in% colnames(metadata[["tx2gene"]]))) {
            message(paste(
                "tx2gene: Renaming `enstxp`, `ensgene`",
                "to `transcriptID`, `geneID`"
            ))
            assert_are_identical(
                x = colnames(metadata[["tx2gene"]]),
                y = c("enstxp", "ensgene")
            )
            colnames(metadata[["tx2gene"]]) <- c("transcriptID", "geneID")
        }
        assertIsTx2gene(metadata[["tx2gene"]])

        # Dead genes: "missing" or "unannotated"
        if ("missingGenes" %in% names(metadata)) {
            message("Dropping missingGenes from metadata")
            metadata[["missingGenes"]] <- NULL
        }
        if ("unannotatedGenes" %in% names(metadata)) {
            message("Dropping unannotatedGenes from metadata")
            metadata[["unannotatedGenes"]] <- NULL
        }

        # yamlFile
        if ("yamlFile" %in% names(metadata)) {
            message("Dropping yamlFile path")
            metadata[["yamlFile"]] <- NULL
        }

        # version
        metadata[["previousVersion"]] <- metadata[["version"]]
        metadata[["version"]] <- packageVersion

        metadata(rse) <- metadata
        rm(metadata)

        # Update assays --------------------------------------------------------
        caller <- metadata(rse)[["caller"]]
        assert_is_a_string(caller)
        assert_is_subset(caller, validCallers)

        level <- metadata(rse)[["level"]]
        assert_is_a_string(level)
        assert_is_subset(level, validLevels)

        # Rename main assay from "raw" to "counts"
        if ("raw" %in% assayNames(rse)) {
            message("Renaming main `raw` assay to `counts`")
            assays(rse)[["counts"]] <- assays(rse)[["raw"]]
            assays(rse)[["raw"]] <- NULL
        }

        # Caller-specific assays
        if (caller %in% tximportCallers) {
            # length (from tximport)
            if (!"length" %in% assayNames(rse)) {
                message("Moving length matrix to assays")
                length <- slot(object, "bcbio")[["tximport"]][["length"]]
                assert_is_matrix(length)
                assays(rse)[["length"]] <- length
            }
            assert_is_subset(tximportAssays, assayNames(rse))
        } else if (caller %in% featureCountsCallers) {
            assert_is_subset(featureCountsAssays, assayNames(rse))
        }

        # Gene-level-specific assays (DESeq2)
        if (level == "genes") {
            # DESeq2 normalized counts
            if (is(assays(rse)[["normalized"]], "DESeqDataSet")) {
                dds <- assays(rse)[["normalized"]]
                assays(rse)[["normalized"]] <- counts(dds, normalized = TRUE)
            } else if (!"normalized" %in% assayNames(rse)) {
                dds <- .regenerateDESeqDataSet(rse)
                assays(rse)[["normalized"]] <- counts(dds, normalized = TRUE)
            }

            # vst
            if (is(assays(rse)[["vst"]], "DESeqTransform")) {
                message("Coercing vst DESeqTransform to matrix in assays")
                assays(rse)[["vst"]] <- assay(assays(rse)[["vst"]])
            }

            # rlog
            if (is(assays(rse)[["rlog"]], "DESeqTransform")) {
                message("Coercing rlog DESeqTransform to matrix in assays")
                assays(rse)[["rlog"]] <- assay(assays(rse)[["rlog"]])
            }
        }

        # Defunct assays
        # tmm
        if ("tmm" %in% assayNames(rse)) {
            message(paste(
                "Dropping tmm from assays.",
                "Now calculated on the fly instead."
            ))
            assays(rse)[["tmm"]] <- NULL
        }

        # Put the required assays first
        assays(rse) <- assays(rse)[unique(c(requiredAssays, assayNames(rse)))]
        assert_is_subset(requiredAssays, assayNames(rse))

        # Column data ----------------------------------------------------------
        colData <- colData(rse)

        # Move metrics from metadata into colData, if necessary
        metrics <- metadata(rse)[["metrics"]]
        if (is.data.frame(metrics)) {
            message("Moving metrics from metadata into colData")

            # Always remove legacy `name` column
            metrics[["name"]] <- NULL

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

            # Only include columns not already present in colData
            setdiff <- setdiff(colnames(metrics), colnames(colData))
            metrics <- metrics[, sort(setdiff), drop = FALSE]

            colData <- cbind(colData, metrics)
            metadata(rse)[["metrics"]] <- NULL
        }

        # Remove legacy `sampleID` and `description` columns, if present
        colData[["sampleID"]] <- NULL
        colData[["description"]] <- NULL

        colData(rse) <- colData

        # Return ---------------------------------------------------------------
        .new.bcbioRNASeq(
            assays = assays(rse),
            # This will handle mismatches
            rowRanges = rowRanges,
            colData = colData(rse),
            metadata = metadata(rse)
        )
    }



#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature("bcbioRNASeq"),
    definition = updateObject.bcbioRNASeq
)
