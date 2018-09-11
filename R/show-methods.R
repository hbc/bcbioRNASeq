#' Show an Object
#'
#' @name show
#' @family S4 Object
#' @author Michael Steinbuagh
#' @importFrom methods show
#' @export
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

        # Include standard `SummarizedExperiment` information.
        return <- c(
            return,
            separatorBar,
            capture.output(show(rse))
        )

        cat(return, sep = "\n")
    }
)



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("DESeqAnalysis"),
    definition = function(object) {
        version <- metadata(object@data)[["version"]]
        contrastNames <- vapply(
            X = object@results,
            FUN = contrastName,
            FUN.VALUE = character(1L)
        )
        return <- c(
            bold(paste(class(object), version)),
            paste(bold("Transform:"), .transformType(object@transform)),
            bold(paste0("Results (", length(contrastNames), "):")),
            paste0("  - ", contrastNames),
            separatorBar,
            capture.output(show(object@data))
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
        up <- slot(object, "degUp")
        down <- slot(object, "degDown")

        contrast <- contrastName(all)
        alpha <- metadata(all)[["alpha"]]
        lfc <- metadata(all)[["lfcThreshold"]]

        summary <- capture.output(summary(all))
        # Remove leading and trailing whitespace.
        summary <- summary[!grepl("^$", summary)]
        # Remove the lines about results documentation.
        summary <- summary[!grepl("\\?results$", summary)]

        # Get back the number of genes that have adjusted P values.
        nPadj <- sum(!is.na(all[["padj"]]))

        # Base mean information.
        nBaseMeanGT1 <- nrow(all[all[["baseMean"]] > 1L, , drop = FALSE])

        return <- c(
            bold(class(object)),
            paste(bold("Contrast:"), contrast),
            paste(bold("Alpha:"), alpha),
            paste(bold("LFC threshold:"), lfc),
            separatorBar,
            paste(nrow(all), "genes total"),
            paste(nBaseMeanGT1, "genes with base mean > 1"),
            paste(nPadj, "genes with adjusted P values"),
            paste(nrow(up), "upregulated genes"),
            paste(nrow(down), "downregulated genes"),
            separatorBar,
            summary
        )

        cat(return, sep = "\n")
    }
)
