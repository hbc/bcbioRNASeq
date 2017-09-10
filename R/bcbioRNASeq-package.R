#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @import methods
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom basejump annotable camel detectOrganism dots fixNA gene2symbol
#'   mdHeader mdList prepareSummarizedExperiment prepareTemplate
#'   readFileByExtension readYAML removeNA snake tx2gene
#' @importFrom BiocGenerics counts design density plotDispEsts plotMA plotPCA
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom data.table rbindlist
#' @importFrom DEGreport degQC degCovariates degPatterns
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix DESeqDataSetFromTximport
#'   DESeqTransform estimateSizeFactors results resultsNames rlog
#'   varianceStabilizingTransformation
#' @importFrom dplyr arrange bind_cols desc distinct group_by left_join mutate
#'   mutate_all mutate_if pull rename ungroup
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom ggplot2 aes_ aes_string annotation_logticks coord_fixed
#'   coord_flip element_blank element_text expand_limits facet_wrap geom_bar
#'   geom_boxplot geom_density geom_hline geom_jitter geom_line geom_point
#'   geom_polygon geom_ribbon geom_smooth geom_text ggplot ggtitle guides labs
#'   scale_color_manual scale_x_continuous scale_x_log10 scale_y_log10 theme
#'   xlab xlim ylab ylim
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom pheatmap pheatmap
#' @importFrom readr read_csv read_delim read_lines read_tsv write_csv
#' @importFrom reshape2 melt
#' @importFrom rlang .data set_names sym syms
#' @importFrom S4Vectors cor head mcols metadata metadata<- na.omit SimpleList
#' @importFrom stats formula setNames
#' @importFrom stringr str_detect str_match str_replace str_replace_all
#' @importFrom SummarizedExperiment assay assays assays<- colData rowData
#'   SummarizedExperiment
#' @importFrom tibble column_to_rownames remove_rownames rownames_to_column
#' @importFrom tidyr expand_
#' @importFrom tximport tximport
#' @importFrom utils capture.output download.file globalVariables packageVersion
#'   read.table
#' @importFrom viridis inferno scale_color_viridis scale_fill_viridis viridis
#' @importFrom vsn meanSdPlot
"_PACKAGE"

globalVariables(".")

# Plot label separator
labelSep <- ": "

projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
metaPriorityCols <- c("sampleID", "sampleName")
perSampleDirs <- c("sailfish", "salmon")

url <- "http://bioinformatics.sph.harvard.edu/bcbioRNASeq"
