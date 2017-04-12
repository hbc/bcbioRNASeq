# General ====
#' @importFrom basejump kables
#' @importFrom knitr kable
#' @importFrom reshape2 melt
#' @importFrom stats cor na.omit
#' @importFrom utils globalVariables head sessionInfo
NULL

#' @importFrom basejump setNamesSnake
set_names_snake <- basejump::setNamesSnake



# RNA-Seq ====
#' @import DESeq2
#' @importFrom biomaRt getBM useEnsembl
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors mcols
#' @importFrom tximport tximport
NULL



# Visualization ====
#' @import ggrepel
#' @importFrom CHBUtils volcano_density_plot
#' @importFrom pheatmap pheatmap
NULL



# GSEA ====
#' @importFrom RDAVIDWebService addList DAVIDWebService getAnnotationSummary
#'   getClusterReport getClusterReportFile getFunctionalAnnotationChart
#'   getFunctionalAnnotationChartFile getFunctionalAnnotationTable
#'   getFunctionalAnnotationTableFile getGeneCategoriesReport getGeneListReport
#'   getGeneListReportFile setTimeOut
NULL



# tidyverse ====
# http://tidyverse.org/
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @import readr
#' @import stringr
#' @import tibble
#' @importFrom tidyr expand_
NULL
