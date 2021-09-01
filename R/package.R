## FIXME Need to import `abort` from basejump, which requires an update.



#' bcbioRNASeq
#'
#' Import and analyze [bcbio](https://bcbio-nextgen.readthedocs.io/) RNA-seq
#' data.
#'
#' @aliases NULL
#' @keywords internal
#'
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom basejump RangedSummarizedExperiment SummarizedExperiment
#'   Tx2Gene
#' @importClassesFrom edgeR DGEList
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom AcidPlots !! acid_coord_flip acid_geom_abline acid_geom_bar
#'   acid_geom_label_repel acid_scale_y_continuous_nopad matchLabels
#'   pretty_breaks plotCountsCorrelation plotCountsCorrelationHeatmap plot_grid
#'   sym
#' @importFrom DESeq2 DESeq DESeqDataSet estimateSizeFactors fpkm rlog
#'   varianceStabilizingTransformation
#' @importFrom basejump DataFrame Rle SimpleList SummarizedExperiment Tx2Gene
#'   alert alertInfo alertSuccess alertWarning as.data.frame assay assay<-
#'   assayNames assayNames<- assays assays<- camelCase capture.output cbind
#'   colData colData<- colSums detectLanes detectOrganism dl do.call droplevels
#'   emptyRanges encode formalsList getListElement h1 h2 h3 head humanize import
#'   importSampleData importTx2Gene interestingGroups interestingGroups<- lapply
#'   makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeSummarizedExperiment mapGenesToRownames matchInterestingGroups mcols
#'   mcols<- metadata metadata<- methodFormals metrics packageName
#'   packageVersion realpath rowData rowRanges rowSums sampleData session_info
#'   showHeader showSlotInfo standardizeCall stripTranscriptVersions txt ul
#'   width
#' @importFrom bcbioBase getGTFFileFromYAML getMetricsFromYAML
#'   getSampleDataFromYAML importDataVersions importProgramVersions projectDir
#'   runDate sampleDirs
#' @importFrom edgeR DGEList calcNormFactors cpm scaleOffset
#' @importFrom ggplot2 aes expand_limits geom_point geom_smooth ggplot ggtitle
#'   guides labs scale_y_continuous theme
#' @importFrom goalie areDisjointSets areIntersectingSets areSetEqual assert
#'   bapply hasDimnames hasLength hasUniqueCols hasValidDimnames isADirectory
#'   isAFile isAURL isAny isCharacter isDirectory isFile isFlag isGGScale isInt
#'   isInRange isNonNegative isNumber isPositive isProportion isString isSubset
#'   validNames validate validateClasses
#' @importFrom methods as as<- is new setAs setClass setMethod setValidity show
#'   slot slot<- validObject .hasSlot
#' @importFrom tximport tximport
"_PACKAGE"
