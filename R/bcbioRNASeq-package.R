#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @name bcbioRNASeq-package
#'
#' @import methods
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#' @importFrom rlang .data abort inform warn
#' @importFrom S4Vectors cor head mcols metadata na.omit SimpleList
#' @importFrom SummarizedExperiment assay assayNames assays colData rowData
#'   rowRanges
#' @importFrom utils globalVariables packageVersion
NULL



globalVariables(".")
packageVersion <- packageVersion("bcbioRNASeq")
metadataPriorityCols <- c("sampleID", "description", "sampleName")
legacyMetricsCols <- c(metadataPriorityCols, "name", "x53Bias")
updateMsg <- "Run `updateObject()` to update your object"
validCallers <- c("salmon", "kallisto", "sailfish")
