#' Melt Count Matrix to Long Format and log10 Transform
#'
#' @rdname meltLog10
#' @name meltLog10
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param normalized Select normalized counts (`TRUE`), raw counts (`FALSE`),
#' or specifically request TMM-normalized counts (`tmm`).
#'
#' @seealso [reshape2::melt()].
#'
#' @return log10 melted [data.frame].
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/rld.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' meltLog10(bcb) %>% glimpse()
#'
#' # DESeqDataSet
#' meltLog10(dds) %>% glimpse()
#'
#' # DESeqTransform
#' meltLog10(rld) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom basejump sanitizeColData
#' @importFrom dplyr left_join
.joinMelt <- function(counts, colData) {
    assert_are_identical(
        colnames(counts),
        as.character(colData[["sampleID"]])
    )
    melted <- .meltLog10(counts)
    colData <- sanitizeColData(colData)
    left_join(melted, colData, by = "sampleID")
}



#' @importFrom dplyr mutate
#' @importFrom magrittr set_colnames
#' @importFrom reshape2 melt
#' @importFrom tibble rownames_to_column
.meltLog10 <- function(counts) {
    assert_is_matrix(counts)
    counts %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        melt(id = 1L) %>%
        set_colnames(c("ensgene", "sampleID", "counts")) %>%
        # Melt operation will define as factor. Let's set this manually later.
        mutate(sampleID = as.character(.data[["sampleID"]])) %>%
        # Remove zero counts
        .[.[["counts"]] > 0L, , drop = FALSE] %>%
        # log10 transform the counts
        mutate(counts = log10(.data[["counts"]])) %>%
        # Arrange by sampleID, to match factor levels
        arrange(.data[["sampleID"]]) %>%
        mutate(sampleID = as.factor(.data[["sampleID"]]))
}



.meltLog10.bcbioRNASeq <- function(  # nolint
    object,
    normalized = TRUE) {
    .joinMelt(
        counts = counts(object, normalized = normalized),
        colData = colData(object)
    )
}



.meltLog10.DESeqDataSet <- function(  # nolint
    object,
    normalized = TRUE) {
    .joinMelt(
        counts = counts(object, normalized = normalized),
        colData = colData(object)
    )
}



.meltLog10.DESeqTransform <- function(object) {  # nolint
    .joinMelt(
        counts = assay(object),
        colData = colData(object)
    )
}



# Methods ======================================================================
#' @rdname meltLog10
#' @export
setMethod(
    "meltLog10",
    signature("bcbioRNASeq"),
    .meltLog10.bcbioRNASeq)



#' @rdname meltLog10
#' @export
setMethod(
    "meltLog10",
    signature("DESeqDataSet"),
    .meltLog10.DESeqDataSet)



#' @rdname meltLog10
#' @export
setMethod(
    "meltLog10",
    signature("DESeqTransform"),
    .meltLog10.DESeqTransform)
