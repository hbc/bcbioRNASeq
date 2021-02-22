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
#'   as.data.frame assay assay<- assayNames assayNames<- assays assays<-
#'   camelCase capture.output cbind colData colData<- colSums detectLanes
#'   detectOrganism do.call droplevels emptyRanges encode formalsList head
#'   humanize import importSampleData importTx2Gene interestingGroups
#'   interestingGroups<- lapply makeGRangesFromEnsembl makeGRangesFromGFF
#'   makeNames makeSummarizedExperiment mapGenesToRownames
#'   matchInterestingGroups mcols mcols<- metadata metadata<- methodFormals
#'   metrics packageName packageVersion realpath rowData rowRanges rowSums
#'   sampleData session_info showSlotInfo standardizeCall
#'   stripTranscriptVersions width
#' @importFrom bcbioBase getGTFFileFromYAML getMetricsFromYAML
#'   getSampleDataFromYAML importDataVersions importProgramVersions projectDir
#'   runDate sampleDirs
#' @importFrom cli cat_line cli_alert cli_alert_info cli_alert_success
#'   cli_alert_warning cli_div cli_dl cli_end cli_h1 cli_h2 cli_h3 cli_text
#'   cli_ul
#' @importFrom edgeR DGEList calcNormFactors cpm scaleOffset
#' @importFrom ggplot2 aes expand_limits geom_point geom_smooth ggplot ggtitle
#'   guides labs scale_y_continuous theme
#' @importFrom goalie areDisjointSets areIntersectingSets areSetEqual assert
#'   bapply containsAURL hasDimnames hasLength hasUniqueCols hasValidDimnames
#'   isADirectory isAFile isAny isCharacter isDirectory isFile isFlag isGGScale
#'   isInt isInRange isNonNegative isNumber isPositive isProportion isString
#'   isSubset validNames validate validateClasses
#' @importFrom methods as as<- is new setAs setClass setMethod setValidity show
#'   slot slot<- validObject .hasSlot
#' @importFrom tximport tximport
"_PACKAGE"



## FIXME Define this in AcidGenerics
#' @importFrom S4Vectors getListElement
NULL
