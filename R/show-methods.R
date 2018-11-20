#' @name show
#' @author Michael Steinbuagh
#' @inherit methods::show
#'
#' @examples
#' data(bcb)
#' show(bcb)
NULL



#' @importFrom methods show
#' @aliases NULL
#' @export
methods::show



.showHeader <- function(object, version = NULL) {
    cat(paste(class(object), version), sep = "\n")
}



# bcbioRNASeq ==================================================================
show.bcbioRNASeq <-  # nolint
    function(object) {
        validObject(object)
        # Metadata.
        m <- metadata(object)
        # Row ranges metadata.
        rrm <- metadata(rowRanges(object))
        .showHeader(object, version = m[["version"]])
        showSlotInfo(list(
            uploadDir = m[["uploadDir"]],
            dates = as.character(c(
                bcbio = m[["runDate"]],
                R = m[["date"]]
            )),
            level = m[["level"]],
            caller = m[["caller"]],
            sampleMetadataFile = m[["sampleMetadataFile"]],
            organism = m[["organism"]],
            gffFile = m[["gffFile"]],
            annotationHub = rrm[["annotationHub"]],
            ensemblRelease = rrm[["release"]],
            genomeBuild = rrm[["build"]],
            interestingGroups = m[["interestingGroups"]]
        ))
        # Extend the standard RangedSummarizedExperiment method.
        rse <- as(object, "RangedSummarizedExperiment")
        cat(capture.output(show(rse)), sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("bcbioRNASeq"),
    definition = show.bcbioRNASeq
)



# DESeqAnalysis ================================================================
show.DESeqAnalysis <-  # nolint
    function(object) {
        validObject(object)
        data <- slot(object, "data")
        transform <- slot(object, "transform")
        .showHeader(
            object = object,
            version = metadata(data)[["version"]]
        )
        contrastNames <- .contrastNames(object)
        showSlotInfo(list(
            transform = .transformType(transform),
            contrastNames = contrastNames
        ))
        cat(capture.output(show(data)), sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("DESeqAnalysis"),
    definition = show.DESeqAnalysis
)



# DESeqResultsTables ===========================================================
show.DESeqResultsTables <-  # nolint
    function(object) {
        validObject(object)
        results <- slot(object, "results")
        deg <- slot(object, "deg")
        metadata <- slot(object, "metadata")

        .showHeader(
            object = object,
            version = metadata[["version"]]
        )

        contrast <- contrastName(results)
        alpha <- metadata(results)[["alpha"]]
        lfcThreshold <- metadata(results)[["lfcThreshold"]]

        list <- list(
            contrast = contrast,
            alpha = alpha,
            lfcThreshold = lfcThreshold
        )

        # Include file paths, if they're stashed (from `export()`).
        if (is.character(metadata[["export"]])) {
            if (isTRUE(metadata[["dropbox"]])) {
                name <- "dropbox"
            } else {
                name <- "dir"
            }
            files <- metadata[["export"]]
            dirname <- unique(dirname(files))
            assert_is_a_string(dirname)
            list[[name]] <- dirname
        }

        showSlotInfo(list)

        # Include DESeqResults summary.
        summary <- capture.output(summary(results)) %>%
            # Remove leading and trailing whitespace.
            .[!grepl("^$", .)] %>%
            # Remove the lines about results documentation.
            .[!grepl("\\?results$", .)]
        cat(summary, sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("DESeqResultsTables"),
    definition = show.DESeqResultsTables
)
