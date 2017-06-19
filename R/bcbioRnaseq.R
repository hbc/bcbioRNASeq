#' bcbioRnaseq
#'
#' Quality control and differential expression for bcbio-nextgen RNA-seq
#' experiments.
#'
#' @keywords internal
#'
#' @import basejump
#' @import Biobase
#' @import DESeq2
#' @import ggplot2
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import rjson
#' @import annotables
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom data.table rbindlist
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom DEGreport degQC degCovariates degPatterns
#' @importFrom ggrepel geom_text_repel
#' @importFrom isomiRs IsomirDataSeqFromFiles
#' @importFrom methods as is new validObject slot slot<-
#' @importFrom pheatmap pheatmap
#' @importFrom stats density
#' @importFrom utils read.table capture.output
#' @importFrom vsn meanSdPlot
"_PACKAGE"

globalVariables(basejump::globals, asNamespace("bcbioRnaseq"), add = TRUE)

fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"

label_sep <- " : "
