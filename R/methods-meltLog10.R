#' Melt Count Matrix to Long Format and log10 Transform
#'
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
.joinMelt <- function(counts, colData) {
    colData <- as.data.frame(colData)
    assert_are_identical(
        colnames(counts),
        as.character(colData[["sampleID"]])
    )
    melted <- .meltLog10(counts)
    assert_is_tbl_df(melted)
    left_join(melted, colData, by = "sampleID")
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
        set_colnames(c("geneID", "sampleID", "counts")) %>%
        .[, c("sampleID", "geneID", "counts")] %>%
        arrange(!!!syms(c("sampleID", "geneID"))) %>%
        group_by(!!sym("sampleID")) %>%
        # Remove zero counts
        .[.[["counts"]] > 0L, , drop = FALSE] %>%
        # log10 transform the counts
        mutate(counts = log10(.data[["counts"]]))
}



# Methods ======================================================================
#' @rdname meltLog10
#' @export
setMethod(
    "meltLog10",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = TRUE
    ) {
        .joinMelt(
            counts = counts(object, normalized = normalized),
            colData = colData(object)
        )
    }
)



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
            colData = colData(object)
        )
    }
)



#' @rdname meltLog10
#' @export
setMethod(
    "meltLog10",
    signature("DESeqTransform"),
    function(object) {
        .joinMelt(
            counts = assay(object),
            colData = colData(object)
        )
    }
)
