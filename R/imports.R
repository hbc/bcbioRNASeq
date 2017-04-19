# General ====
#' @importFrom knitr asis_output kable opts_knit
#' @importFrom reshape2 melt
#' @importFrom stats cor density na.omit
#' @importFrom utils globalVariables head sessionInfo
NULL



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
#' @import grid
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom pheatmap pheatmap
NULL



# tidyverse ====
# http://tidyverse.org/
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @import readr
#' @import stringr
#' @import tibble
#' @importFrom rlang .data UQ
#' @importFrom tidyr expand_
NULL
