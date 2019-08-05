#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#' @importClassesFrom basejump Tx2Gene
#' @importClassesFrom edgeR DGEList
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom BiocGenerics cbind colSums density design lapply sapply
#'   updateObject width
#' @importFrom DEGreport degCovariates
#' @importFrom DESeq2 DESeq DESeqDataSet DESeqDataSetFromMatrix
#'   DESeqDataSetFromTximport DESeqResults DESeqTransform estimateDispersions
#'   estimateSizeFactors fpkm priorInfo results resultsNames rlog
#'   varianceStabilizingTransformation
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom S4Vectors DataFrame Rle as.data.frame complete.cases decode
#'   getListElement head mcols mcols<- metadata metadata<- na.omit
#' @importFrom SummarizedExperiment SummarizedExperiment assay assay<-
#'   assayNames assayNames<- assays assays<- colData colData<- rowData rowRanges
#' @importFrom acidplots acid_geom_abline acid_geom_label acid_geom_label_repel
#'   acid_coord_flip acid_geom_bar acid_scale_y_continuous_nopad
#'   plotCorrelationHeatmap plotCountsCorrelation plotCountsCorrelationHeatmap
#'   plotFeaturesDetected plotHeatmap plotPCA
#' @importFrom basejump Gene2Symbol Tx2Gene camelCase coerceS4ToList
#'   convertGenesToSymbols detectLanes detectOrganism emptyRanges formalsList
#'   import initDir interestingGroups interestingGroups<- lanePattern
#'   makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeSummarizedExperiment mapGenesToRownames markdownHeader markdownList
#'   markdownPlotlist matchArgsToDoCall matchInterestingGroups
#'   matchesGene2Symbol meltCounts methodFormals metrics organism
#'   prepareTemplate printString realpath relevel removeNA sampleData
#'   sampleData<- separator showSlotInfo standardizeCall stripTranscriptVersions
#'   uniteInterestingGroups
#' @importFrom bcbioBase copyToDropbox getGTFFileFromYAML getMetricsFromYAML
#'   getSampleDataFromYAML metadataBlacklist projectDir projectDirPattern
#'   readDataVersions readProgramVersions readSampleData readTx2Gene runDate
#'   sampleDirs
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_cols bind_rows desc everything filter group_by
#'   left_join mutate mutate_all mutate_if pull rename row_number select
#'   select_if starts_with ungroup
#' @importFrom edgeR DGEList calcNormFactors cpm scaleOffset
#' @importFrom ggplot2 aes annotation_logticks coord_fixed coord_flip
#'   element_blank element_text expand_limits facet_wrap geom_bar geom_boxplot
#'   geom_density geom_hline geom_jitter geom_point geom_polygon geom_ribbon
#'   geom_smooth geom_violin geom_vline ggplot ggtitle guides labs
#'   position_jitterdodge scale_color_manual scale_x_continuous scale_x_discrete
#'   scale_y_continuous stat_summary theme xlab ylim
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom goalie areDisjointSets areIntersectingSets assert bapply
#'   containsAURL hasDimnames hasLength hasRownames hasUniqueCols
#'   hasValidDimnames isADirectory isAFile isAny isCharacter isDirectory isFile
#'   isFlag isGGScale isInt isInRange isNonEmpty isNonNegative isNumber
#'   isPositive isProportion isString isSubset validNames validate
#'   validateClasses
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom methods as as<- is new setAs setClass setMethod setValidity show
#'   slot slot<- validObject .hasSlot
#' @importFrom readr read_csv read_tsv write_csv
#' @importFrom rlang !! !!! := UQ sym syms
#' @importFrom scales pretty_breaks
#' @importFrom sessioninfo session_info
#' @importFrom stringr str_match str_trunc
#' @importFrom tibble as_tibble column_to_rownames remove_rownames
#'   rownames_to_column tibble
#' @importFrom tximport tximport
#' @importFrom utils capture.output globalVariables packageVersion
#' @importFrom vsn meanSdPlot
NULL
