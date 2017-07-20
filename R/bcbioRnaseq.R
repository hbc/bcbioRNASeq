#' bcbioRnaseq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @import basejump BiocGenerics DESeq2 SummarizedExperiment S4Vectors rjson
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
#' @importFrom isomiRs IsomirDataSeqFromFiles
#' @importFrom methods as is new slot slot<- validObject
#' @importFrom pheatmap pheatmap
#' @importFrom stats formula
#' @importFrom tximport tximport
#' @importFrom utils read.table capture.output
#' @importFrom vsn meanSdPlot
"_PACKAGE"

globalVariables(".")

# Quality control plot colors
qc_fail_color <- "red"
qc_pass_color <- "green"
qc_warn_color <- "orange"
qc_line_alpha <- 0.75
qc_line_size <- 2L

# Plot label separator
label_sep <- ": "

project_dir_pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
meta_priority_cols <- c("sample_id", "sample_name")
per_sample_dirs <- c("sailfish", "salmon")
