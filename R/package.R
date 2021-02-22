#' bcbioRNASeq
#'
#' Import and analyze [bcbio](https://bcbio-nextgen.readthedocs.io/) RNA-seq
#' data.
#'
#' @aliases NULL
#' @keywords internal
#'
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'   SummarizedExperiment
#' @importClassesFrom basejump Tx2Gene
#' @importClassesFrom edgeR DGEList
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom AcidPlots acid_coord_flip acid_geom_abline acid_geom_bar
#'   acid_geom_label_repel acid_scale_y_continuous_nopad matchLabels
#'   plotCountsCorrelation plotCountsCorrelationHeatmap
#' @importFrom BiocGenerics colSums rowSums updateObject
#' @importFrom DESeq2 DESeq DESeqDataSet estimateSizeFactors fpkm rlog
#'   varianceStabilizingTransformation
#' @importFrom S4Vectors DataFrame Rle SimpleList as.data.frame cbind do.call
#'   getListElement head lapply mcols mcols<- metadata metadata<- sapply width
#' @importFrom SummarizedExperiment SummarizedExperiment assay assay<-
#'   assayNames assayNames<- assays assays<- colData colData<- rowData rowRanges
#' @importFrom basejump Tx2Gene camelCase detectLanes detectOrganism droplevels
#'   emptyRanges encode formalsList humanize import importSampleData
#'   importTx2Gene interestingGroups interestingGroups<- makeGRangesFromEnsembl
#'   makeGRangesFromGFF makeNames makeSummarizedExperiment mapGenesToRownames
#'   matchInterestingGroups methodFormals metrics realpath sampleData
#'   showSlotInfo standardizeCall stripTranscriptVersions
#' @importFrom bcbioBase getGTFFileFromYAML getMetricsFromYAML
#'   getSampleDataFromYAML importDataVersions importProgramVersions projectDir
#'   runDate sampleDirs
#' @importFrom cli cat_line cli_alert cli_alert_info cli_alert_success
#'   cli_alert_warning cli_div cli_dl cli_end cli_h1 cli_h2 cli_h3 cli_text
#'   cli_ul
#' @importFrom cowplot plot_grid
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
#' @importFrom rlang !! sym
#' @importFrom scales pretty_breaks
#' @importFrom sessioninfo session_info
#' @importFrom tximport tximport
#' @importFrom utils capture.output packageName packageVersion
"_PACKAGE"
