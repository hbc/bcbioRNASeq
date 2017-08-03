#' bcbioRnaseq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @import basejump BiocGenerics DESeq2 rjson SummarizedExperiment S4Vectors
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom data.table rbindlist
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom DEGreport degQC degCovariates degPatterns
#' @importFrom ggplot2 aes_ coord_fixed coord_flip element_blank element_text
#'   expand_limits facet_wrap geom_bar geom_boxplot geom_density geom_hline
#'   geom_jitter geom_line geom_point geom_polygon geom_ribbon geom_smooth
#'   ggplot ggtitle guides labs scale_x_continuous scale_y_log10 theme xlab xlim
#'   ylab ylim scale_x_log10 scale_color_manual geom_text aes_string
#'   scale_x_log10 annotation_logticks
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom isomiRs IsomirDataSeqFromFiles
#' @importFrom methods as is new slot slot<- validObject
#' @importFrom pheatmap pheatmap
#' @importFrom stats formula
#' @importFrom tximport tximport
#' @importFrom utils read.table capture.output
#' @importFrom vsn meanSdPlot
"_PACKAGE"

globalVariables(".")

packageURL <- "http://bioinformatics.sph.harvard.edu/bcbioRnaseq"

# Quality control plot colors
qcFailColor <- "red"
qcPassColor <- "green"
qcWarnColor <- "orange"
qcLineAlpha <- 0.75
qcLineSize <- 2L

# Plot label separator
labelSep <- ": "

projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
metaPriorityCols <- c("sampleID", "sampleName")
perSampleDirs <- c("sailfish", "salmon")
