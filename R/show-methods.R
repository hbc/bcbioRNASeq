#' Show an object
#' @name show
#' @author Michael Steinbuagh
#' @inherit methods::show
#' @note Updated 2019-08-07.
#' @examples
#' data(bcb)
#' show(bcb)
NULL



## Updated 2019-07-23.
.showHeader <- function(object, version = NULL) {
    cat(paste(class(object), version), sep = "\n")
}



## Updated 2019-07-23.
`show,bcbioRNASeq` <-  # nolint
    function(object) {
        validObject(object)
        ## Metadata.
        m <- metadata(object)
        ## Row ranges metadata.
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
        ## Extend the standard RangedSummarizedExperiment method.
        rse <- as(object, "RangedSummarizedExperiment")
        cat(capture.output(show(rse)), sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("bcbioRNASeq"),
    definition = `show,bcbioRNASeq`
)
