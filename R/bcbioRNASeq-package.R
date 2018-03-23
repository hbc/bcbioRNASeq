#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @name bcbioRNASeq-package
#'
#' @import S4Vectors methods
#'
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#'
#' @importFrom SummarizedExperiment assay assayNames assays colData rowData
#'   rowRanges
#' @importFrom basejump camel
#' @importFrom dplyr filter left_join pull
#' @importFrom ggplot2 aes_string element_text expand_limits geom_jitter ggplot
#'   ggtitle labs theme
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom rlang !!! !! abort inform sym syms warn
#' @importFrom tibble rownames_to_column
#' @importFrom utils globalVariables packageVersion
NULL



globalVariables(".")
packageVersion <- packageVersion("bcbioRNASeq")
lanePattern <- "_L(\\d{3})"
metadataPriorityCols <- c("sampleID", "description", "sampleName")
legacyMetricsCols <- c(metadataPriorityCols, "name", "x53Bias")
updateMsg <- "Run `updateObject()` to update your object"
validCallers <- c("salmon", "kallisto", "sailfish")
requiredAssays <- c("raw", "tpm", "length")
