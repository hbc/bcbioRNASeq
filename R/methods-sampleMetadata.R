#' Sample Metadata
#'
#' Returns the sample metadata.
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#'
#' @importFrom bcbioBase sampleMetadata
#'
#' @inheritParams general
#'
#' @return [data.frame].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' sampleMetadata(bcb) %>% glimpse()
#'
#' # DESeqDataSet
#' dds <- bcbio(bcb, "DESeqDataSet")
#' sampleMetadata(dds) %>% glimpse()
#'
#' # DESeqTransform
#' rld <- assays(bcb)[["rlog"]]
#' sampleMetadata(rld) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr mutate_all
#' @importFrom magrittr set_rownames
.sampleMetadata <- function(object, ...) {
    data <- colData(object)
    assert_is_all_of(data, "DataFrame")
    assert_is_subset("sampleID", colnames(data))
    data %>%
        as.data.frame() %>%
        mutate_all(as.factor) %>%
        set_rownames(.[["sampleID"]])
}



# Methods ======================================================================
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
