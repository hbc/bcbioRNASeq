#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#' @importClassesFrom basejump Tx2Gene
#' @importClassesFrom edgeR DGEList
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom BiocGenerics colSums rowSums updateObject
#' @importFrom DESeq2 DESeq DESeqDataSet estimateSizeFactors fpkm rlog
#'   varianceStabilizingTransformation
#' @importFrom S4Vectors DataFrame Rle SimpleList as.data.frame cbind do.call
#'   getListElement head lapply mcols mcols<- metadata metadata<- sapply width
#' @importFrom SummarizedExperiment SummarizedExperiment assay assay<-
#'   assayNames assayNames<- assays assays<- colData colData<- rowData rowRanges
#' @importFrom acidplots acid_coord_flip acid_geom_abline acid_geom_bar
#'   acid_geom_label_repel acid_scale_y_continuous_nopad matchLabels
#'   plotCountsCorrelation plotCountsCorrelationHeatmap
#' @importFrom basejump Tx2Gene camelCase detectLanes detectOrganism droplevels
#'   emptyRanges encode formalsList humanize import importSampleData
#'   importTx2Gene interestingGroups interestingGroups<- makeGRangesFromEnsembl
#'   makeGRangesFromGFF makeNames makeSummarizedExperiment mapGenesToRownames
#'   matchArgsToDoCall matchInterestingGroups methodFormals metrics realpath
#'   sampleData showSlotInfo standardizeCall stripTranscriptVersions
#' @importFrom bcbioBase getGTFFileFromYAML getMetricsFromYAML
#'   getSampleDataFromYAML importDataVersions importProgramVersions projectDir
#'   runDate sampleDirs
#' @importFrom cowplot plot_grid
#' @importFrom edgeR DGEList calcNormFactors cpm scaleOffset
#' @importFrom ggplot2 aes expand_limits geom_point geom_smooth ggplot ggtitle
#'   guides labs scale_y_continuous theme
#' @importFrom goalie areDisjointSets areIntersectingSets areSetEqual assert
#'   bapply containsAURL hasDimnames hasLength hasUniqueCols hasValidDimnames
#'   isADirectory isAFile isAny isCharacter isDirectory isFile isFlag isGGScale
#'   isInt isInRange isNonEmpty isNonNegative isNumber isPositive isProportion
#'   isString isSubset validNames validate validateClasses
#' @importFrom methods as as<- is new setAs setClass setMethod setValidity show
#'   slot slot<- validObject .hasSlot
#' @importFrom rlang !! sym
#' @importFrom scales pretty_breaks
#' @importFrom sessioninfo session_info
#' @importFrom tximport tximport
#' @importFrom utils capture.output globalVariables packageVersion
NULL
