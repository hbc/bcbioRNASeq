#' bcbioRNASeq
#'
#' Quality control and differential expression for
#' [bcbio](http://bcbio-nextgen.readthedocs.io) RNA-seq experiments.
#'
#' @aliases NULL
#' @keywords internal
#'
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#' @importClassesFrom basejump Tx2Gene
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom BiocGenerics cbind colSums density design lapply sapply width
#' @importFrom DESeq2 DESeq DESeqDataSet DESeqDataSetFromMatrix
#'   DESeqDataSetFromTximport DESeqResults DESeqTransform estimateDispersions
#'   estimateSizeFactors fpkm priorInfo results resultsNames rlog
#'   varianceStabilizingTransformation
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom S4Vectors DataFrame Rle as.data.frame complete.cases decode head
#'   mcols mcols<- metadata metadata<- na.omit
#' @importFrom SummarizedExperiment SummarizedExperiment assay assay<-
#'   assayNames assayNames<- assays assays<- colData colData<- rowData rowRanges
#' @importFrom basejump Gene2Symbol Tx2Gene basejump_geom_abline
#'   basejump_geom_label basejump_geom_label_repel camel coerceS4ToList
#'   convertGenesToSymbols detectLanes detectOrganism emptyRanges formalsList
#'   import initDir interestingGroups interestingGroups<- lanePattern
#'   makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeSummarizedExperiment mapGenesToRownames markdownHeader markdownList
#'   markdownPlotlist matchArgsToDoCall matchInterestingGroups
#'   matchesGene2Symbol meltCounts methodFormals metrics organism
#'   plotGenesDetected plotHeatmap prepareTemplate printString realpath
#'   relevelColData relevelRowRanges removeNA sampleData sampleData<-
#'   sanitizeRowData sanitizeSampleData separator showSlotInfo snake
#'   standardizeCall stripTranscriptVersions uniteInterestingGroups
#' @importFrom bcbioBase copyToDropbox getGTFFileFromYAML getMetricsFromYAML
#'   getSampleDataFromYAML metadataBlacklist projectDir projectDirPattern
#'   readDataVersions readProgramVersions readSampleData readTx2Gene runDate
#'   sampleDirs
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_cols desc everything filter group_by left_join
#'   mutate mutate_all mutate_if pull rename row_number select select_if
#'   starts_with
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom ggplot2 aes annotation_logticks coord_fixed coord_flip
#'   element_blank element_text expand_limits facet_wrap geom_bar geom_boxplot
#'   geom_density geom_hline geom_jitter geom_point geom_polygon geom_ribbon
#'   geom_smooth geom_violin geom_vline ggplot ggtitle guides labs
#'   position_jitterdodge scale_color_manual scale_x_continuous
#'   scale_y_continuous stat_summary theme xlab ylim
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom goalie areDisjointSets areIntersectingSets assert hasDimnames
#'   hasLength hasRownames hasUniqueCols hasValidDimnames isADirectory isAFile
#'   isAny isCharacter isDirectory isFile isFlag isGGScale isInt isInRange
#'   isNonEmpty isNonNegative isNumber isPositive isProportion isString isSubset
#'   validNames validate validateClasses
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom matrixStats colMedians
#' @importFrom methods .hasSlot as as<- is new setAs show slot slot<-
#'   validObject
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
"_PACKAGE"
