#' \code{bcbioRnaseq} package
#'
#' Quality control and differential expression for bcbio-nextgen RNA-seq
#' experiments.
#'
#' See the README on \href{https://github.com/roryk/bcbioRnaseq}{GitHub}.
#'
#' @docType package
#' @name bcbioRnaseq
NULL



# Globals ====
globalVariables(".")
pass_color <- "green"
warn_color <- "orange"
fail_color <- "red"



# Imports ====
# General ----
#' @importFrom knitr asis_output kable opts_knit
#' @importFrom reshape2 melt
#' @importFrom stats cor density na.omit
#' @importFrom utils globalVariables head sessionInfo
NULL

# RNA-Seq ----
#' @import DESeq2
#' @importFrom biomaRt getBM useEnsembl
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors mcols
#' @importFrom tximport tximport
NULL

# Visualization ----
#' @import ggrepel
#' @importFrom pheatmap pheatmap
NULL

# tidyverse ----
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



# Rexports ====
#' @keywords internal
#' @usage NULL
#' @export
magrittr::`%>%`
