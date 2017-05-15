#' bcbioRnaseq
#'
#' Quality control and differential expression for bcbio-nextgen RNA-seq
#' experiments.
#'
#' @import basejump
#' @import DESeq2
#' @import ggplot2
#' @import methods
#' @import SummarizedExperiment
#' @import S4Vectors
#' @importFrom cowplot ggdraw draw_plot
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom DEGreport degQC
#' @importFrom ggrepel geom_text_repel
#' @importFrom isomiRs IsomirDataSeqFromFiles
#' @importFrom pheatmap pheatmap
#' @importFrom stats density
#' @importFrom tximport tximport
#' @importFrom utils read.table
#' @importFrom vsn meanSdPlot
"_PACKAGE"

globalVariables(basejump::globals,
                asNamespace("bcbioRnaseq"),
                add = TRUE)

fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"
