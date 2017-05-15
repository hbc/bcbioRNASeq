#' bcbioRnaseq.
#'
#' Quality control and differential expression for bcbio-nextgen RNA-seq
#' experiments.
#'
#' @import basejump
#' @import DESeq2
#' @import ggplot2
#' @importFrom biomaRt getBM listMarts useEnsembl
#' @importFrom cowplot ggdraw draw_plot
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom DEGreport degQC
#' @importFrom ggrepel geom_text_repel
#' @importFrom isomiRs IsomirDataSeqFromFiles
#' @importFrom methods as is new validObject
#' @importFrom pheatmap pheatmap
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors mcols
#' @importFrom tximport tximport
#' @importFrom vsn meanSdPlot
#' @importFrom yaml yaml.load_file
"_PACKAGE"

globalVariables(basejump::globals,
                asNamespace("bcbioRnaseq"),
                add = TRUE)

fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"
