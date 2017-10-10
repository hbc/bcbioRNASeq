#' Sample Metadata
#'
#' Returns the sample metadata.
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame].
#'
#' @examples
#' data(bcb, dds, rld)
#'
#' # bcbioRNASeq
#' sampleMetadata(bcb)
#'
#' # DESeqDataSet
#' sampleMetadata(dds)
#'
#' # DESeqTransform
#' sampleMetadata(dds)
NULL



# Constructors ====
.sampleMetadata <- function(object, ...) {
    object %>%
        colData() %>%
        as.data.frame()
}



# Methods ====
#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata",
    signature("bcbioRNASeqANY"),
    .sampleMetadata)



#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata",
    signature("DESeqANY"),
    .sampleMetadata)
