#' Update object
#'
#' @name updateObject
#' @author Michael Steinbaugh
#' @note Updated 2022-04-25.
#'
#' @details
#' Update old objects created by the bcbioRNASeq package. The session
#' information metadata is preserved from the time when the bcbio data was
#' originally loaded into R.
#'
#' @section Legacy `bcbioRNADataSet` class:
#'
#' Support for `bcbioRNADataSet` objects was dropped in v0.2.0 of the package.
#' If you need to load one of these objects, please install an older release.
#'
#' @section Legacy `bcbioRnaseq` package:
#'
#' The previous `bcbioRnaseq` package (note case) must be reinstalled to load
#' objects from versions <= 0.0.20. We changed the name of the package to
#' `bcbioRNASeq` starting in v0.0.21.
#'
#' @inheritParams AcidRoxygen::params
#' @param rowRanges `GRanges` or `NULL`.
#' Row annotations. Since we converted to `RangedSummarizedExperiment` in
#' v0.2.0, this option had to be added to enable updating of newly required
#' `rowRanges` slot. Objects that are >= v0.2 don't require this argument and
#' it can be left `NULL`.
#'
#' @return Modified object.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' object <- bcb
#' object <- updateObject(object)
#' validObject(object)
#'
#' ## Example that depends on remote file.
#' ## > object <- import(file.path(
#' ## >     bcbioRNASeqTestsURL,
#' ## >     "bcbioRNASeq_0.1.4.rds"
#' ## > ))
#' ## > object <- updateObject(object)
#' ## > validObject(object)
NULL



## Updated 2022-04-25.
`updateObject,bcbioRNASeq` <- # nolint
    function(object,
             rowRanges = NULL,
             ...,
             verbose = FALSE) {
        assert(isFlag(verbose))
        metadata <- metadata(object)
        ## Renamed 'version' to 'packageVersion' on 2021-03-16.
        version <- metadata[["packageVersion"]]
        if (is.null(version)) {
            version <- metadata[["version"]]
        }
        assert(is(version, "package_version"))
        if (isTRUE(verbose)) {
            h1("Update object")
            alert(sprintf(
                fmt = paste(
                    "Updating {.var %s} object from version {.val %s}",
                    "to {.val %s}."
                ),
                "bcbioRNASeq",
                as.character(version),
                as.character(.pkgVersion)
            ))
        }
        ## Legacy slots --------------------------------------------------------
        if (isTRUE(verbose)) {
            h2("Legacy slots")
        }
        ## NAMES
        if (!is.null(slot(object, "NAMES"))) {
            if (isTRUE(verbose)) {
                alertInfo(sprintf(
                    "{.var %s} slot must be set to {.val %s}.",
                    "NAMES", "NULL"
                ))
            }
            slot(object, "NAMES") <- NULL
        }
        ## elementMetadata
        if (ncol(slot(object, "elementMetadata")) != 0L) {
            if (isTRUE(verbose)) {
                alertInfo(sprintf(
                    fmt = paste(
                        "{.var %s} slot must contain a",
                        "zero-column {.var %s}."
                    ),
                    "elementMetadata", "DataFrame"
                ))
            }
            slot(object, "elementMetadata") <-
                as(matrix(nrow = nrow(object), ncol = 0L), "DataFrame")
        }
        ## rowRanges
        if (!.hasSlot(object, "rowRanges")) {
            if (is.null(rowRanges)) {
                if (isTRUE(verbose)) {
                    alertWarning(sprintf(
                        "Slotting empty {.fun %s}.", "rowRanges"
                    ))
                }
                assays <- slot(object, "assays")
                ## Extract assay matrix from ShallowSimpleListAssays object.
                if (packageVersion("SummarizedExperiment") >= "1.15") {
                    ## Bioconductor 3.10+.
                    assay <- getListElement(x = assays, i = 1L)
                } else {
                    ## Legacy method.
                    assay <- assays[[1L]]
                }
                rownames <- rownames(assay)
                assert(isCharacter(rownames))
                rowRanges <- emptyRanges(names = rownames)
                rowData <- slot(object, "elementMetadata")
                mcols(rowRanges) <- rowData
            }
            assert(isAny(rowRanges, c("GenomicRanges", "GenomicRangesList")))
            slot(object, "rowRanges") <- rowRanges
        } else if (
            .hasSlot(object, "rowRanges") &&
                !is.null(rowRanges)
        ) {
            abort(sprintf(
                fmt = paste(
                    "Object already contains {.fun %s}.",
                    "Don't attempt to slot new ones with {.arg %s} argument.",
                    sep = "\n"
                ),
                "rowRanges", "rowRanges"
            ))
        }
        ## Check for legacy bcbio slot.
        if (.hasSlot(object, "bcbio")) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Dropping legacy {.var %s} slot.", "bcbio"
                ))
            }
        }
        ## Metadata ------------------------------------------------------------
        if (isTRUE(verbose)) {
            h2("Metadata")
        }
        ## bcbioLog
        if (is.null(metadata[["bcbioLog"]])) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "bcbioLog", "empty character"
                ))
            }
            metadata[["bcbioLog"]] <- character()
        }
        ## bcbioCommandsLog
        if (is.null(metadata[["bcbioCommandsLog"]])) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "bcbioCommands", "empty character"
                ))
            }
            metadata[["bcbioCommandsLog"]] <- character()
        }
        ## call
        if (!"call" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Stashing empty {.var %s}.", "call"
                ))
            }
            metadata[["call"]] <- call(name = "bcbioRNASeq")
        }
        ## caller
        if (!"caller" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "caller", "salmon"
                ))
            }
            metadata[["caller"]] <- "salmon"
        }
        ## dataVersions
        dataVersions <- metadata[["dataVersions"]]
        if (is(dataVersions, "data.frame")) {
            metadata[["dataVersions"]] <- as(dataVersions, "DataFrame")
        }
        ## design
        if ("design" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Dropping legacy {.var %s}.", "design"
                ))
            }
            metadata[["design"]] <- NULL
        }
        ## ensemblRelease
        if ("ensemblVersion" %in% names(metadata)) {
            ## Renamed in v0.2.0.
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "ensemblVersion", "ensemblRelease"
                ))
            }
            names(metadata)[
                names(metadata) == "ensemblVersion"
            ] <- "ensemblRelease"
        }
        if (!is.integer(metadata[["ensemblRelease"]])) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "ensemblRelease", "integer"
                ))
            }
            metadata[["ensemblRelease"]] <-
                as.integer(metadata[["ensemblRelease"]])
        }
        ## genomeBuild
        if (!is.character(metadata[["genomeBuild"]])) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "genomeBuild", "empty character"
                ))
            }
            metadata[["genomeBuild"]] <- character()
        }
        ## gffFile
        if ("gtfFile" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "gtfFile", "gffFile"
                ))
            }
            names(metadata)[
                names(metadata) == "gtfFile"
            ] <- "gffFile"
        }
        if (!"gffFile" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "gffFile", "empty character"
                ))
            }
            metadata[["gffFile"]] <- character()
        }
        ## gtf
        if ("gtf" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Dropping {.val %s} in {.var %s}.",
                    "stashed GTF", "gtf"
                ))
            }
            metadata[["gtf"]] <- NULL
        }
        ## lanes
        if (!is.integer(metadata[["lanes"]])) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as integer.", "lanes"
                ))
            }
            metadata[["lanes"]] <- as.integer(metadata[["lanes"]])
        }
        ## level
        if (!"level" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "level", "genes"
                ))
            }
            metadata[["level"]] <- "genes"
        }
        ## packageVersion
        metadata[["previousVersion"]] <- metadata[["version"]]
        metadata[["packageVersion"]] <- .pkgVersion
        metadata[["version"]] <- NULL
        ## programVersions
        if (!"programVersions" %in% names(metadata) &&
            "programs" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "programs", "programVersions"
                ))
            }
            names(metadata)[
                names(metadata) == "programs"
            ] <- "programVersions"
        }
        programVersions <- metadata[["programVersions"]]
        if (is(programVersions, "data.frame")) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Coercing {.var %s} to {.var %s}.",
                    "programVersions", "DataFrame"
                ))
            }
            metadata[["programVersions"]] <- as(programVersions, "DataFrame")
        }
        ## sampleMetadataFile
        if (!is.character(metadata[["sampleMetadataFile"]])) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "sampleMetadataFile", "empty character"
                ))
            }
            metadata[["sampleMetadataFile"]] <- character()
        }
        ## sessionInfo
        ## Support for legacy devtoolsSessionInfo stash.
        ## Previously, we stashed both devtools* and utils* variants.
        if ("devtoolsSessionInfo" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Simplifying stashed {.var %s}.", "sessionInfo"
                ))
            }
            names(metadata)[
                names(metadata) == "devtoolsSessionInfo"
            ] <- "sessionInfo"
            metadata[["utilsSessionInfo"]] <- NULL
        }
        ## template
        if ("template" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alert(sprintf("Dropping legacy {.var %s}.", "template"))
            }
            metadata[["template"]] <- NULL
        }
        ## Dead genes: "missing" or "unannotated".
        if ("missingGenes" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Dropping {.var %s} from {.fun %s}.",
                    "missingGenes", "metadata"
                ))
            }
            metadata[["missingGenes"]] <- NULL
        }
        if ("unannotatedGenes" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alertWarning(sprintf(
                    "Dropping {.var %s} from {.fun %s}.",
                    "unannotatedGenes", "metadata"
                ))
            }
            metadata[["unannotatedGenes"]] <- NULL
        }
        ## yamlFile
        if ("yamlFile" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Dropping {.var %s} file path.", "yamlFile"
                ))
            }
            metadata[["yamlFile"]] <- NULL
        }
        ## tximport-specific
        if (isSubset(metadata[["caller"]], .tximportCallers)) {
            if (isTRUE(verbose)) {
                h3("tximport")
            }
            ## countsFromAbundance
            if (!"countsFromAbundance" %in% names(metadata)) {
                countsFromAbundance <- "lengthScaledTPM"
                if (isTRUE(verbose)) {
                    alertWarning(sprintf(
                        "Setting {.var %s} as {.val %s}.",
                        "countsFromAbundance", countsFromAbundance
                    ))
                }
                metadata[["countsFromAbundance"]] <- countsFromAbundance
            }
            ## tx2gene
            tx2gene <- metadata[["tx2gene"]]
            if (!is(tx2gene, "Tx2Gene")) {
                if (isTRUE(verbose)) {
                    alert(sprintf(
                        "Coercing {.var %s} to {.var %s} class object.",
                        "tx2gene", "Tx2Gene"
                    ))
                }
                assert(is.data.frame(tx2gene))
                tx2gene <- as(tx2gene, "DataFrame")
                colnames(tx2gene) <- c("txId", "geneId")
                rownames(tx2gene) <- NULL
                metadata[["tx2gene"]] <- Tx2Gene(tx2gene)
            }
        }
        ## Filter NULL metadata.
        metadata <- Filter(f = Negate(is.null), x = metadata)
        ## Assays --------------------------------------------------------------
        if (isTRUE(verbose)) {
            h2("Assays")
        }
        assays <- assays(object)
        ## These variables are needed for assay handling.
        caller <- metadata[["caller"]]
        assert(
            isString(caller),
            isSubset(caller, .callers)
        )
        level <- metadata[["level"]]
        assert(
            isString(level),
            isSubset(level, .levels)
        )
        ## Ensure raw counts are always named "counts".
        if ("raw" %in% names(assays)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} assay to {.var %s}.",
                    "raw", "counts"
                ))
            }
            names(assays)[names(assays) == "raw"] <- "counts"
        }
        ## Rename average transcript length matrix.
        if ("length" %in% names(assays)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} assay to {.var %s}.",
                    "length", "avgTxLength"
                ))
            }
            names(assays)[names(assays) == "length"] <- "avgTxLength"
        }
        ## Drop legacy TMM counts.
        if ("tmm" %in% names(assays)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    fmt = paste(
                        "Dropping {.var %s} from {.fun %s}.",
                        "Calculating on the fly instead."
                    ),
                    "tmm", "assays"
                ))
            }
            assays[["tmm"]] <- NULL
        }
        ## Gene-level-specific assays (DESeq2). Handle legacy objects where size
        ## factor normalized counts aren't stashed as a matrix. Note that these
        ## will always be length scaled.
        if (level == "genes") {
            ## DESeq2 normalized counts.
            if (is(assays[["normalized"]], "DESeqDataSet")) {
                if (isTRUE(verbose)) {
                    alert(sprintf(
                        fmt = paste(
                            "Coercing {.var %s} assay from",
                            "{.val %s} to {.val %s}."
                        ),
                        "normalized", "DESeqDataSet", "matrix"
                    ))
                }
                dds <- assays[["normalized"]]
                assays[["normalized"]] <- counts(dds, normalized = TRUE)
            }
            ## Variance-stabilizing transformation.
            if (is(assays[["vst"]], "DESeqTransform")) {
                if (isTRUE(verbose)) {
                    alert(sprintf(
                        "Coercing {.var %s} assay from {.val %s} to {.val %s}.",
                        "vst", "DESeqTransform", "matrix"
                    ))
                }
                assays[["vst"]] <- assay(assays[["vst"]])
            }
            ## Regularized log.
            if (is(assays[["rlog"]], "DESeqTransform")) {
                if (isTRUE(verbose)) {
                    alert(sprintf(
                        "Coercing {.var %s} assay from {.val %s} to {.val %s}.",
                        "rlog", "DESeqTransform", "matrix"
                    ))
                }
                assays[["rlog"]] <- assay(assays[["rlog"]])
            }
        }
        ## Always put the required assays first.
        assert(isSubset(.assays, names(assays)))
        assays <- assays[unique(c(.assays, names(assays)))]
        ## Check for required caller-specific assays.
        if (caller %in% .tximportCallers) {
            assert(isSubset(.tximportAssays, names(assays)))
        } else if (caller %in% .featureCountsCallers) {
            assert(isSubset(.featureCountsAssays, names(assays)))
        }
        ## Filter NULL assays.
        assays <- Filter(f = Negate(is.null), x = assays)
        ## Row ranges ----------------------------------------------------------
        if (isTRUE(verbose)) {
            h2("Row ranges")
        }
        rowRanges <- rowRanges(object)
        if (hasLength(colnames(mcols(rowRanges)))) {
            colnames(mcols(rowRanges)) <-
                camelCase(colnames(mcols(rowRanges)), strict = TRUE)
            if (any(grepl(
                pattern = "^transcript",
                x = colnames(mcols(rowRanges))
            ))) {
                colnames(mcols(rowRanges)) <- gsub(
                    pattern = "^transcript",
                    replacement = "tx",
                    x = colnames(mcols(rowRanges))
                )
            }
        }
        ## rowRangesMetadata
        if ("rowRangesMetadata" %in% names(metadata)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Moving {.val %s} into {.var %s} metadata.",
                    "rowRangesMetadata", "rowRanges"
                ))
            }
            metadata(rowRanges)[["ensembldb"]] <-
                metadata[["rowRangesMetadata"]]
            metadata[["rowRangesMetadata"]] <- NULL
        }
        ## biotype
        if ("biotype" %in% colnames(mcols(rowRanges))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "biotype", "geneBiotype"
                ))
            }
            mcols(rowRanges)[["geneBiotype"]] <-
                as.factor(mcols(rowRanges)[["biotype"]])
            mcols(rowRanges)[["biotype"]] <- NULL
        }
        ## broadClass
        if (
            "broadClass" %in% colnames(mcols(rowRanges)) &&
                is.character(mcols(rowRanges)[["broadClass"]])
        ) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} to {.val %s}.",
                    "broadClass", "factor"
                ))
            }
            mcols(rowRanges)[["broadClass"]] <-
                as.factor(mcols(rowRanges)[["broadClass"]])
        }
        ## ensgene
        if ("ensgene" %in% colnames(mcols(rowRanges))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "ensgene", "geneId"
                ))
            }
            mcols(rowRanges)[["geneId"]] <-
                as.character(mcols(rowRanges)[["ensgene"]])
            mcols(rowRanges)[["ensgene"]] <- NULL
        }
        ## symbol
        if ("symbol" %in% colnames(mcols(rowRanges))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "symbol", "geneName"
                ))
            }
            mcols(rowRanges)[["geneName"]] <-
                as.factor(mcols(rowRanges)[["symbol"]])
            mcols(rowRanges)[["symbol"]] <- NULL
        }
        ## Use run-length encoding (Rle) for metadata columns.
        ## Recommended method as of v0.3.0 update.
        rowRanges <- encode(rowRanges)
        ## Column data ---------------------------------------------------------
        if (isTRUE(verbose)) {
            h2("Column data")
        }
        colData <- colData(object)
        colnames(colData) <- camelCase(colnames(colData), strict = TRUE)
        ## Move metrics from metadata into colData, if necessary.
        metrics <- metadata[["metrics"]]
        if (!is.null(metrics)) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Moving {.var %s} from {.fun %s} into {.fun %s}.",
                    "metrics", "metadata", "colData"
                ))
            }
            assert(is.data.frame(metrics))
            ## Always remove legacy name column.
            metrics[["name"]] <- NULL
            ## Rename 5'3' bias.
            if ("x53Bias" %in% colnames(metrics)) {
                if (isTRUE(verbose)) {
                    alert(sprintf(
                        "Renaming {.var %s} to {.var %s}.",
                        "x53Bias", "x5x3Bias"
                    ))
                }
                metrics[["x5x3Bias"]] <- metrics[["x53Bias"]]
                metrics[["x53Bias"]] <- NULL
            }
            ## Rename rRNA rate.
            if (!"rrnaRate" %in% colnames(metrics)) {
                col <- grep(
                    pattern = "rrnarate",
                    x = colnames(metrics),
                    ignore.case = TRUE,
                    value = TRUE
                )
                assert(isString(col))
                if (isTRUE(verbose)) {
                    alert(sprintf(
                        "Renaming {.var %s} to {.var %s}.",
                        col, "rrnaRate"
                    ))
                }
                metrics[["rrnaRate"]] <- metrics[[col]]
                metrics[[col]] <- NULL
            }
            ## Only include columns not already present in colData.
            setdiff <- setdiff(colnames(metrics), colnames(colData))
            metrics <- metrics[, sort(setdiff), drop = FALSE]
            colData <- cbind(colData, metrics)
            metadata[["metrics"]] <- NULL
        }
        ## Remove legacy sample identifier and description columns, if present.
        colData[["description"]] <- NULL
        colData[["sampleID"]] <- NULL
        colData[["sampleId"]] <- NULL
        ## Return --------------------------------------------------------------
        se <- SummarizedExperiment(
            assays = assays,
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
        bcb <- new(Class = "bcbioRNASeq", se)
        validObject(bcb)
        if (isTRUE(verbose)) {
            alertSuccess(sprintf(
                "Update of {.var %s} object was successful.",
                "bcbioRNASeq"
            ))
        }
        bcb
    }



#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature(object = "bcbioRNASeq"),
    definition = `updateObject,bcbioRNASeq`
)
