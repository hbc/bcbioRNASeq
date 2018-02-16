#' Sample Metadata
#'
#' Returns the sample metadata.
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#'
#' @importFrom bcbioBase sampleMetadata sampleMetadata<-
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
#'
#' # Assignment support
NULL



# Constructors =================================================================
.sampleMetadata <- function(object, ...) {
    as.data.frame(colData(object))
}



`.sampleMetadata<-` <- function(object, ..., value) {
    colData(object) <- value
    object
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



# Assignment Methods ===========================================================
#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata<-",
    signature(
        object = "bcbioRNASeq",
        value = "data.frame"
    ),
    `.sampleMetadata<-`)




#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata<-",
    signature(
        object = "bcbioRNASeq",
        value = "DataFrame"
    ),
    `.sampleMetadata<-`)



#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata<-",
    signature(
        object = "DESeqDataSet",
        value = "data.frame"
    ),
    `.sampleMetadata<-`)



#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata<-",
    signature(
        object = "DESeqDataSet",
        value = "DataFrame"
    ),
    `.sampleMetadata<-`)



#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata<-",
    signature(
        object = "DESeqTransform",
        value = "data.frame"
    ),
    `.sampleMetadata<-`)



#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata<-",
    signature(
        object = "DESeqTransform",
        value = "DataFrame"
    ),
    `.sampleMetadata<-`)
