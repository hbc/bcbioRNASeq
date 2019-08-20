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
#' @importFrom DESeq2 DESeq DESeqDataSet DESeqDataSetFromMatrix
#'   DESeqDataSetFromTximport DESeqResults DESeqTransform estimateDispersions
#'   estimateSizeFactors fpkm priorInfo results resultsNames rlog
#'   varianceStabilizingTransformation
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom S4Vectors DataFrame Rle SimpleList as.data.frame complete.cases
#'   decode getListElement head mcols mcols<- metadata metadata<- na.omit
#' @importFrom SummarizedExperiment SummarizedExperiment assay assay<-
#'   assayNames assayNames<- assays assays<- colData colData<- rowData rowRanges
#' @importFrom acidplots acid_geom_abline acid_geom_label acid_geom_label_repel
#'   acid_coord_flip acid_geom_bar acid_scale_y_continuous_nopad
#'   plotCorrelationHeatmap plotCountsCorrelation plotCountsCorrelationHeatmap
#'   plotFeaturesDetected plotHeatmap plotPCA
#' @importFrom basejump Gene2Symbol Tx2Gene camelCase coerceS4ToList
#'   convertGenesToSymbols detectLanes detectOrganism droplevels emptyRanges
#'   encode formalsList import initDir interestingGroups interestingGroups<-
#'   lanePattern makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeSummarizedExperiment mapGenesToRownames markdownHeader markdownList
#'   markdownPlotlist matchArgsToDoCall matchInterestingGroups
#'   matchesGene2Symbol meltCounts methodFormals metrics organism
#'   prepareTemplate printString readSampleData readTx2Gene realpath removeNA
#'   sampleData sampleData<- separator showSlotInfo standardizeCall
#'   stripTranscriptVersions uniteInterestingGroups
#' @importFrom bcbioBase getGTFFileFromYAML getMetricsFromYAML
#'   getSampleDataFromYAML projectDir readDataVersions readProgramVersions
#'   runDate sampleDirs
#' @importFrom cowplot draw_plot ggdraw plot_grid
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
#' @importFrom methods as as<- is new setAs setClass setMethod setValidity show
#'   slot slot<- validObject .hasSlot
#' @importFrom rlang !! sym
#' @importFrom scales pretty_breaks
#' @importFrom sessioninfo session_info
#' @importFrom stringr str_match str_trunc
#' @importFrom tximport tximport
#' @importFrom utils capture.output globalVariables packageVersion
NULL
