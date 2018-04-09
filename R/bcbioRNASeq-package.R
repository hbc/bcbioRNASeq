#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @name bcbioRNASeq-package
#' @keywords internal
#'
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#' @importFrom BiocGenerics colSums density design
#' @importFrom DEGreport degCovariates
#' @importFrom DESeq2 DESeq DESeqDataSetFromTximport DESeqTransform
#'   estimateSizeFactors results rlog varianceStabilizingTransformation
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom S4Vectors as.data.frame head mcols metadata na.omit
#' @importFrom SummarizedExperiment assay assayNames assays assays<- colData
#'   colData<- rowData rowRanges SummarizedExperiment
#' @importFrom basejump camel convertGenesToSymbols detectOrganism fixNA
#'   initializeDirectory makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeTx2geneFromGFF markdownHeader markdownList markdownPlotlist readYAML
#'   sanitizeRowData sanitizeSampleData snake
#' @importFrom bcbioBase copyToDropbox flatFiles gene2symbol interestingGroups
#'   plotHeatmap prepareSummarizedExperiment prepareTemplate readDataVersions
#'   readLog readProgramVersions readSampleData readTx2gene sampleData
#'   sampleDirs sampleYAMLMetadata sampleYAMLMetrics uniteInterestingGroups
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_cols desc filter group_by left_join mutate
#'   mutate_all mutate_if pull select_if
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @importFrom ggplot2 aes_ aes_string annotation_logticks coord_fixed
#'   coord_flip element_blank element_text expand_limits facet_wrap geom_bar
#'   geom_boxplot geom_density geom_hline geom_jitter geom_point geom_polygon
#'   geom_ribbon geom_smooth ggplot ggtitle guides labs position_jitterdodge
#'   scale_color_manual scale_fill_manual scale_x_continuous scale_x_log10
#'   stat_summary theme xlab ylim
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom methods .hasSlot as is new show slot validObject
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
