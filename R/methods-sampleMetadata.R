#' Sample Metadata
#'
#' Returns the sample metadata.
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#'
#' @importFrom basejump sampleMetadata
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
    signature("bcbioRNASeq"),
    .sampleMetadata)



#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata",
    signature("DESeqDataSet"),
    .sampleMetadata)



#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata",
    signature("DESeqTransform"),
    .sampleMetadata)
