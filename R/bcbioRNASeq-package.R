#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @name bcbioRNASeq-package
#'
#' @import S4Vectors methods
#'
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#'
#' @importFrom BiocGenerics density design
#' @importFrom DEGreport degCovariates
#' @importFrom DESeq2 DESeq DESeqDataSetFromTximport DESeqTransform
#'   estimateSizeFactors results rlog varianceStabilizingTransformation
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom SummarizedExperiment assay assayNames assays colData rowData
#'   rowRanges SummarizedExperiment
#' @importFrom basejump camel convertGenesToSymbols detectOrganism ensembl fixNA
#'   gene2symbol initializeDirectory makeNames markdownHeader markdownList
#'   markdownPlotlist readYAML sanitizeRowData sanitizeSampleData snake
#' @importFrom bcbioBase copyToDropbox interestingGroups
#'   prepareSummarizedExperiment prepareTemplate readDataVersions readLogFile
#'   readProgramVersions readSampleMetadataFile sampleDirs sampleYAMLMetadata
#'   sampleYAMLMetrics uniteInterestingGroups
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_cols desc filter group_by left_join mutate
#'   mutate_all mutate_if pull select_if
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom ggplot2 aes_ aes_string annotation_logticks coord_fixed
#'   coord_flip element_blank element_text expand_limits geom_bar geom_boxplot
#'   geom_density geom_hline geom_jitter geom_point geom_polygon geom_ribbon
#'   geom_smooth ggplot ggtitle guides labs position_jitterdodge
#'   scale_color_manual scale_fill_manual scale_x_continuous scale_x_log10
#'   stat_summary theme xlab ylim
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom parallel mcmapply
#' @importFrom readr read_csv read_tsv write_csv
#' @importFrom reshape2 melt
#' @importFrom rlang !! !!! abort inform sym syms warn
#' @importFrom stringr str_match str_trunc
#' @importFrom tibble as_tibble column_to_rownames glimpse remove_rownames
#'   rownames_to_column tibble
#' @importFrom tximport tximport
#' @importFrom utils capture.output globalVariables packageVersion
#' @importFrom viridis inferno scale_color_viridis scale_fill_viridis viridis
#' @importFrom vsn meanSdPlot
NULL



globalVariables(".")
packageVersion <- packageVersion("bcbioRNASeq")
lanePattern <- "_L(\\d{3})"
metadataPriorityCols <- c("sampleID", "description", "sampleName")
legacyMetricsCols <- c(metadataPriorityCols, "name", "x53Bias")
updateMsg <- "Run `updateObject()` to update your object"
validCallers <- c("salmon", "kallisto", "sailfish")
requiredAssays <- c("raw", "tpm", "length")
