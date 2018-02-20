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
#' @importFrom dplyr left_join
.joinMelt <- function(counts, metadata) {
    assert_are_identical(
        colnames(counts),
        as.character(metadata[["sampleID"]])
    )
    melted <- .meltLog10(counts)
    assert_is_tbl_df(melted)
    assert_is_data.frame(metadata)
    left_join(melted, metadata, by = "sampleID")
}



#' @importFrom dplyr group_by mutate
#' @importFrom magrittr set_colnames
#' @importFrom reshape2 melt
#' @importFrom rlang !! !!! sym syms
#' @importFrom tibble as_tibble rownames_to_column
.meltLog10 <- function(counts) {
    assert_is_matrix(counts)
    counts %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        melt(id = 1L) %>%
        as_tibble() %>%
        set_colnames(c("ensgene", "sampleID", "counts")) %>%
        .[, c("sampleID", "ensgene", "counts")] %>%
        arrange(!!!syms(c("sampleID", "ensgene"))) %>%
        group_by(!!sym("sampleID")) %>%
        # Remove zero counts
        .[.[["counts"]] > 0L, , drop = FALSE] %>%
        # log10 transform the counts
        mutate(counts = log10(.data[["counts"]]))
}



.meltLog10.bcbioRNASeq <- function(  # nolint
    object,
    normalized = TRUE) {
    .joinMelt(
        counts = counts(object, normalized = normalized),
        metadata = sampleMetadata(object)
    )
}



.meltLog10.DESeqDataSet <- function(  # nolint
    object,
    normalized = TRUE) {
    .joinMelt(
        counts = counts(object, normalized = normalized),
        metadata = sampleMetadata(object)
    )
}



.meltLog10.DESeqTransform <- function(object) {  # nolint
    .joinMelt(
        counts = assay(object),
        metadata = sampleMetadata(object)
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
