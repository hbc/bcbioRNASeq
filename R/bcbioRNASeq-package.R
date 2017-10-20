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
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom utils globalVariables
NULL

utils::globalVariables(".")

projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
perSampleDirs <- c("sailfish", "salmon")
