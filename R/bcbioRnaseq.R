#' bcbioRnaseq
#'
#' Quality control and differential expression for bcbio-nextgen RNA-seq
#' experiments. Consult the
#' \href{http://bioinformatics.sph.harvard.edu/bcbioRnaseq}{package website} for
#' additional information.
#'
#' @docType package
#' @name bcbioRnaseq
NULL



# Globals ====
globalVariables(".")
fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"



# Imports ====
## General ----
#' @importFrom knitr asis_output kable opts_chunk opts_knit
#' @importFrom reshape2 melt
#' @importFrom stats cor density na.omit
#' @importFrom utils download.file globalVariables head read.table sessionInfo
#' @importFrom yaml yaml.load_file
NULL

## RNA-Seq ----
#' @import DESeq2
#' @importFrom biomaRt getBM listMarts useEnsembl
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom DEGreport degQC
#' @importFrom RDAVIDWebService DAVIDWebService addList getAnnotationSummary
#'     getClusterReport getClusterReportFile getFunctionalAnnotationChart
#'     getFunctionalAnnotationChartFile getFunctionalAnnotationTable
#'     getFunctionalAnnotationTableFile getGeneCategoriesReport
#'     getGeneListReport getGeneListReportFile setTimeOut
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors mcols
#' @importFrom tximport tximport
#' @importFrom vsn meanSdPlot
NULL

## small RNA-Seq ----
#' @importFrom isomiRs IsomirDataSeqFromFiles


## Visualization ----
#' @import ggrepel
#' @importFrom pheatmap pheatmap
#' @importFrom cowplot ggdraw draw_plot
NULL

## tidyverse ----
## http://tidyverse.org/
#' @import dplyr
#' @import ggplot2
#' @import readr
#' @import readxl
#' @import stringr
#' @importFrom magrittr %>% set_names set_rownames
#' @importFrom rlang !!! !! .data sym syms UQ
#' @importFrom tibble as_tibble glimpse is_tibble remove_rownames rownames_to_column
#' @importFrom tidyr expand_
NULL



# Re-exports ====
## RNA-seq ----
#' @usage NULL
#' @export
DESeq2::counts

#' @usage NULL
#' @export
DESeq2::design

#' @usage NULL
#' @export
vsn::meanSdPlot

## Quasi-quotation ----
#' @usage NULL
#' @export
rlang::`!!!`

#' @usage NULL
#' @export
rlang::`!!`

#' @usage NULL
#' @export
rlang::sym

#' @usage NULL
#' @export
rlang::syms

## Knit utilities ----
#' @usage NULL
#' @export
magrittr::`%>%`

#' @usage NULL
#' @export
tibble::glimpse

#' @usage NULL
#' @export
knitr::kable

#' @usage NULL
#' @export
knitr::opts_chunk

#' @usage NULL
#' @export
magrittr::set_rownames
