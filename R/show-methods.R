#' @name show
#' @author Michael Steinbuagh
#' @importFrom methods show
#' @inherit methods::show
#' @export
#'
#' @examples
#' data(bcb_small)
#' show(bcb_small)
NULL



.showHeader <- function(object, version = NULL) {
    cat(c(
        bold(paste(class(object), version)),
        italic("http://bioinformatics.sph.harvard.edu/bcbioRNASeq"),
        "citation(\"bcbioRNASeq\")"
    ), sep = "\n")
}


.show.bcbioRNASeq <-  # nolint
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
        cat(capture.output(show(as(object, "RangedSummarizedExperiment"))))
    }



.show.DESeqAnalysis <-  # nolint
    function(object) {
        validObject(object)
        .showHeader(
            object = object,
            version = metadata(object@data)[["version"]]
        )
        contrastNames <- vapply(
            X = object@results,
            FUN = contrastName,
            FUN.VALUE = character(1L)
        )
        showSlotInfo(list(
            transform = .transformType(object@transform),
            contrastNames = contrastNames
        ))
        cat(capture.output(show(object@data)))
    }



.show.DESeqResultsTables <-  # nolint
    function(object) {
        validObject(object)
        .showHeader(object)

        all <- slot(object, "all")
        up <- slot(object, "degUp")
        down <- slot(object, "degDown")

        contrast <- contrastName(all)
        alpha <- metadata(all)[["alpha"]]
        lfcThreshold <- metadata(all)[["lfcThreshold"]]

        list <- list(
            contrast = contrast,
            alpha = alpha,
            lfcThreshold = lfcThreshold
        )

        # Include file paths, if they're stashed (from `export()`).
        if (length(object@dropboxFiles) > 0L) {
            name <- "dropbox export"
            paths <- object@dropboxFiles
        } else if (length(object@localFiles) > 0L) {
            name <- "export"
            paths <- object@localFiles
        } else {
            paths <- character()
        }
        if (length(paths) > 0L) {
            dirname <- unique(dirname(paths))
            assert_is_a_string(dirname)
            list[[name]] <- dirname
        }

        showSlotInfo(list)

        # Include DESeqResults summary.
        summary <- capture.output(summary(all)) %>%
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
    signature = signature("bcbioRNASeq"),
    definition = .show.bcbioRNASeq
)



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("DESeqAnalysis"),
    definition = .show.DESeqAnalysis
)



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("DESeqResultsTables"),
    definition = .show.DESeqResultsTables
)
