#' Show an Object
#'
#' @name show
#' @family S4 Object
#' @author Michael Steinbuagh
#'
#' @inherit methods::show
#'
#' @examples
#' show(bcb_small)
NULL



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("bcbioRNASeq"),
    definition = function(object) {
        validObject(object)

        # Extend the RangedSummarizedExperiment method
        rse <- as(object, "RangedSummarizedExperiment")

        return <- c(
            bold(paste(class(object), metadata(object)[["version"]])),
            "http://bioinformatics.sph.harvard.edu/bcbioRNASeq",
            "citation(\"bcbioRNASeq\")",
            separatorBar,
            capture.output(show(rse)),
            separatorBar,
            paste(
                bold("Upload Dir:"),
                deparse(metadata(object)[["uploadDir"]])
            ),
            paste(
                bold("Upload Date:"),
                metadata(object)[["runDate"]]
            ),
            paste(
                bold("R Load Date:"),
                metadata(object)[["date"]]
            ),
            paste(
                bold("Level:"),
                deparse(metadata(object)[["level"]])
            ),
            paste(
                bold("Caller:"),
                deparse(metadata(object)[["caller"]])
            ),
            paste(
                bold("Organism:"),
                deparse(metadata(object)[["organism"]])
            ),
            paste(
                bold("Interesting Groups:"),
                deparse(metadata(object)[["interestingGroups"]])
            )
        )

        # sampleMetadataFile
        sampleMetadataFile <- metadata(object)[["sampleMetadataFile"]]
        if (length(sampleMetadataFile)) {
            return <- c(
                return,
                paste(bold("Metadata File:"), deparse(sampleMetadataFile))
            )
        }

        # Gene annotations
        m <- metadata(object)[["rowRangesMetadata"]]
        if (is.data.frame(m) && length(m)) {
            annotationHub <-
                m[m[["name"]] == "id", "value", drop = TRUE]
            ensemblRelease <-
                m[m[["name"]] == "ensembl_version", "value", drop = TRUE]
            genomeBuild <-
                m[m[["name"]] == "genome_build", "value", drop = TRUE]
            return <- c(
                return,
                paste(bold("AnnotationHub:"), deparse(annotationHub)),
                paste(bold("Ensembl Release:"), deparse(ensemblRelease)),
                paste(bold("Genome Build:"), deparse(genomeBuild))
            )
        }

        # GFF File
        gffFile <- metadata(object)[["gffFile"]]
        if (length(gffFile)) {
            return <- c(
                return,
                paste(bold("GFF File:"), deparse(gffFile))
            )
        }

        cat(return, sep = "\n")
    }
)



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("DESeqAnalysis"),
    definition = function(object) {
        version <- metadata(object@DESeqDataSet)[["version"]]
        contrastNames <- vapply(
            X = object@DESeqResults,
            FUN = contrastName,
            FUN.VALUE = character(1L)
        )
        return <- c(
            bold(paste(class(object), version)),
            printString(slotNames(object)),
            separatorBar,
            capture.output(show(object@DESeqDataSet)),
            separatorBar,
            paste(bold("Transform:"), .transformType(object@DESeqTransform)),
            bold(paste0("Results (", length(contrastNames), "):")),
            paste0("  - ", contrastNames)
        )
        cat(return, sep = "\n")
    }
)



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("DESeqResultsTables"),
    definition = function(object) {
        validObject(object)

        all <- slot(object, "all")
        up <- slot(object, "upregulated")
        down <- slot(object, "downregulated")

        contrast <- contrastName(all)
        alpha <- metadata(all)[["alpha"]]
        lfc <- metadata(all)[["lfcThreshold"]]

        summary <- capture.output(summary(all))
        # Remove leading and trailing whitespace.
        summary <- summary[!grepl("^$", summary)]

        # Get back the number of genes that have adjusted P values.
        nPadj <- sum(!is.na(all[["padj"]]))

        return <- c(
            bold(class(object)),
            paste(bold("Contrast:"), contrast),
            paste(bold("Alpha:"), alpha),
            paste(bold("LFC threshold:"), lfc),
            separatorBar,
            paste(nrow(object@all), "genes total"),
            paste(nPadj, "genes with adjusted P values"),
            paste(nrow(up), "upregulated genes"),
            paste(nrow(down), "downregulated genes"),
            separatorBar,
            summary
        )

        cat(return, sep = "\n")
    }
)
