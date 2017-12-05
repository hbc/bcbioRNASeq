#' Melt Count Matrix to Long Format and log10 Transform
#'
#' @rdname meltLog10
#' @name meltLog10
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param normalized Select normalized counts (`TRUE`), raw counts (`FALSE`),
#' or specifically request TMM-normalized counts (`tmm`).
#'
#' @seealso [reshape2::melt()].
#'
#' @return log10 melted [data.frame].
#'
#' @examples
#' # bcbioRNASeq
#' bcb <- examples[["bcb"]]
#' meltLog10(bcb) %>% glimpse()
#'
#' # DESeqDataSet
#' dds <- examples[["dds"]]
#' meltLog10(dds) %>% glimpse()
#'
#' # DESeqTransform
#' rld <- examples[["rld"]]
#' meltLog10(rld) %>% glimpse()
NULL



# Constructors ====
#' @importFrom dplyr left_join
.joinMelt <- function(counts, metadata) {
    if (!identical(
            colnames(counts),
            as.character(metadata[["sampleID"]])
        )) {
        stop("Sample name mismatch between counts and metadata",
             call. = FALSE)
    }
    melted <- .meltLog10(counts)
    metadata <- as.data.frame(metadata) %>%
        mutate_if(is.factor, droplevels)
    left_join(melted, metadata, by = "sampleID")
}



#' @importFrom dplyr mutate
#' @importFrom reshape2 melt
#' @importFrom stats setNames
#' @importFrom tibble rownames_to_column
.meltLog10 <- function(counts) {
    counts %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        melt(id = 1) %>%
        setNames(c("ensgene", "sampleID", "counts")) %>%
        # Melt operation will define as factor. Let's set this manually later.
        mutate(sampleID = as.character(.data[["sampleID"]])) %>%
        # Remove zero counts
        .[.[["counts"]] > 0, , drop = FALSE] %>%
        # log10 transform the counts
        mutate(counts = log10(.data[["counts"]])) %>%
        # Arrange by sampleID, to match factor levels
        arrange(.data[["sampleID"]]) %>%
        mutate(sampleID = as.factor(.data[["sampleID"]]))
}



# Methods ====
#' @rdname meltLog10
#' @export
setMethod(
    "meltLog10",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = TRUE) {
        .joinMelt(
            counts = counts(object, normalized = normalized),
            metadata = colData(object))
    })



#' @rdname meltLog10
#' @export
setMethod(
    "meltLog10",
    signature("DESeqDataSet"),
    function(
        object,
        normalized = TRUE) {
        .joinMelt(
            counts = counts(object, normalized = normalized),
            metadata = colData(object))
    })



#' @rdname meltLog10
#' @export
setMethod(
    "meltLog10",
    signature("DESeqTransform"),
    function(object) {
        .joinMelt(
            counts = assay(object),
            metadata = colData(object))
    })
