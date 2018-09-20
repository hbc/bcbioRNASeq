setClassUnion(
    name = "missingOrNULL",
    members = c("missing", "NULL")
)



# bcbioRNASeq ==================================================================
#' `bcbioRNASeq` Class
#'
#' `bcbioRNASeq` is an S4 class that extends `RangedSummarizedExperiment`, and
#' is designed to store a [bcbio](https://bcbio-nextgen.readthedocs.org) RNA-seq
#' analysis.
#'
#' @note `bcbioRNASeq` extended `SummarizedExperiment` prior to v0.2.0, where we
#'   migrated to `RangedSummarizedExperiment`.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh, Lorena Pantano
#' @export
setClass(
    Class = "bcbioRNASeq",
    contains = "RangedSummarizedExperiment"
)



# DESeqAnalysis ================================================================
#' `DESeqAnalysis` Class
#'
#' Class containing all elements generated during differential expression
#' analysis with DESeq2. This class is essentially a `list` with validity checks
#' to ensure `DESeqTransform` and `DESeqResults` correspond to the
#' `DESeqDataSet`.
#'
#' @section DESeqDataSet:
#'
#' We recommend generating the `DESeqDataSet` by coercion from `bcbioRNASeq`
#' object using `as(dds, "bcbioRNASeq")`. Don't use the [DESeq2::DESeqDataSet()]
#' or [DESeq2::DESeqDataSetFromMatrix()] constructors to generate the
#' `DESeqDataSet` object.
#'
#' @section DESeqResults:
#'
#' Don't modify any of the `DESeqResults` objects manually. This includes
#' rearranging the rows or dropping genes without adjusted P values. We'll take
#' care of this automatically in supported methods.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @slot data `DESeqDataSet`.
#' @slot transform `DESeqTransform`.
#' @slot results `list`. One or more unshrunken `DESeqResults`.
#' @slot lfcShrink `list`. One or more shrunken `DESeqResults`.
setClass(
    Class = "DESeqAnalysis",
    slots = c(
        data = "DESeqDataSet",
        transform = "DESeqTransform",
        results = "list",
        lfcShrink = "list"
    ),
    prototype = list(
        lfcShrink = list()
    )
)



# DESeqResultsTables ===========================================================
#' `DESeqResultsTables` Class
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @slot all `DESeqResults`. Original unmodified `DESeqResults`. Should contain
#'   all genes, including those with `NA` adjusted *P* values.
#' @slot deg `DataFrame`. Subset containing genes that pass adjusted *P* value
#'   and log2 fold change cutoffs.
#' @slot degUp `DataFrame`. Directional subset containing only upregulated
#'   genes.
#' @slot degDown `DataFrame`. Directional subset containing only downregulated
#'   genes.
#' @slot localFiles `list`. Local file paths.
#' @slot dropboxFiles `list`. Dropbox file paths.
setClass(
    Class = "DESeqResultsTables",
    slots = c(
        all = "DESeqResults",
        deg = "DataFrame",
        degUp = "DataFrame",
        degDown = "DataFrame",
        localFiles = "character",
        dropboxFiles = "list"
    ),
    # Consider setting an initialize method instead.
    prototype = list(
        localFiles = character(),
        dropboxFiles = list()
    )
)
