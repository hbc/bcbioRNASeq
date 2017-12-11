#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @rdname bcbioRNASeq-package
#' @name bcbioRNASeq-package
#'
#' @import methods SummarizedExperiment
#'
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
NULL

#' @importFrom utils globalVariables
globalVariables(".")

#' @importFrom utils packageVersion
packageVersion <- packageVersion("bcbioRNASeq")

projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
metadataPriorityCols <- c("sampleID", "description", "sampleName")
perSampleDirs <- c("sailfish", "salmon")
