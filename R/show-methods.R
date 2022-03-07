#' Show an object
#'
#' @name show
#' @author Michael Steinbaugh
#' @note Updated 2022-03-07.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return Console output.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' show(bcb)
NULL



## Updated 2021-02-22.
`show,bcbioRNASeq` <-  # nolint
    function(object) {
        validObject(object)
        showHeader(object)
        ## Metadata.
        m <- metadata(object)
        ## Row ranges metadata.
        rrm <- metadata(rowRanges(object))
        list <- list(
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
        )
        if (.isFastMode(object)) {
            list[["fast"]] <- TRUE
        }
        showSlotInfo(list)
        ## Extend the standard RangedSummarizedExperiment method.
        rse <- as(object, "RangedSummarizedExperiment")
        cat(capture.output(show(rse)), sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "bcbioRNASeq"),
    definition = `show,bcbioRNASeq`
)
