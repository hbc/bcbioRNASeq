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



# Methods ======================================================================
#' @rdname show
#' @export
setMethod(
    "show",
    signature("bcbioRNASeq"),
    function(object) {
        validObject(object)

        # Extend the RangedSummarizedExperiment show method
        rse <- as(object, "RangedSummarizedExperiment")
        return <- capture.output(show(rse))

        return[[1L]] <- paste(
            "class:",
            class(object),
            metadata(object)[["version"]]
        )

        return <- c(
            return,
            bcbioBase::separatorBar,
            paste(
                "Upload Dir:",
                deparse(metadata(object)[["uploadDir"]])
            ),
            paste(
                "Upload Date:",
                metadata(object)[["runDate"]]
            ),
            paste(
                "R Load Date:",
                metadata(object)[["date"]]
            ),
            paste(
                "Level:",
                deparse(metadata(object)[["level"]])
            ),
            paste(
                "Caller:",
                deparse(metadata(object)[["caller"]])
            ),
            paste(
                "Organism:",
                deparse(metadata(object)[["organism"]])
            )
        )

        # sampleMetadataFile
        sampleMetadataFile <- metadata(object)[["sampleMetadataFile"]]
        if (length(sampleMetadataFile)) {
            return <- c(
                return,
                paste("Metadata File:", deparse(sampleMetadataFile))
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
                paste("AnnotationHub:", deparse(annotationHub)),
                paste("Ensembl Release:", deparse(ensemblRelease)),
                paste("Genome Build:", deparse(genomeBuild))
            )
        }

        # GFF File
        gffFile <- metadata(object)[["gffFile"]]
        if (length(gffFile)) {
            return <- c(
                return,
                paste("GFF File:", deparse(gffFile))
            )
        }

        return <- c(
            return,
            paste("Installed Version:", packageVersion("bcbioRNASeq"))
        )

        cat(return, sep = "\n")
    }
)
