## S4 classes ==================================================================

#' @importClassesFrom AcidGenomes Tx2Gene
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom IRanges DataFrameList
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' SummarizedExperiment
#' @importClassesFrom edgeR DGEList
NULL



## S4 generics and methods =====================================================

#' @importFrom AcidGenerics as.DESeqDataSet as.DESeqTransform as.DGEList
#' camelCase droplevels2 encode humanize interestingGroups interestingGroups<-
#' makeNames mapGenesToRownames markdownTables metrics plot5Prime3PrimeBias
#' plotCorrelationHeatmap plotCounts plotCountsCorrelation
#' plotCountsCorrelationHeatmap plotCountsPerBiotype plotCountsPerBroadClass
#' plotCountsPerFeature plotDEGHeatmap plotExonicMappingRate plotGenderMarkers
#' plotGeneSaturation plotIntronicMappingRate plotMappedReads plotMappingRate
#' plotPCACovariates plotMeanSD plotPseudoVsAlignedCounts plotQC
#' plotRRNAMappingRate plotTotalReads relativeLogExpression sampleData
#' showHeader slotAlignedCounts stripTranscriptVersions tmm topTables
#' @importFrom AcidGenomes Tx2Gene
#' @importFrom BiocGenerics as.data.frame colSums counts counts<- do.call
#' estimateSizeFactors plotDispEsts plotPCA rowSums updateObject width
#' @importFrom S4Vectors Rle cbind droplevels getListElement head lapply mcols
#' mcols<- metadata metadata<-
#' @importFrom SummarizedExperiment assay assay<- assayNames assayNames<- assays
#' assays<- colData colData<- rowData rowRanges
#' @importFrom methods coerce show
#' @importFrom pipette export import
#'
#' @importMethodsFrom AcidBase showHeader
#' @importMethodsFrom AcidExperiment droplevels2 humanize interestingGroups
#' interestingGroups<- mapGenesToRownames metrics sampleData
#' stripTranscriptVersions
#' @importMethodsFrom AcidGenomes stripTranscriptVersions
#' @importMethodsFrom AcidPlots plotCorrelationHeatmap plotCountsCorrelation
#' plotCountsCorrelationHeatmap plotCountsPerBiotype plotCountsPerBroadClass
#' plotCountsPerFeature plotGenderMarkers plotPCA
#' @importMethodsFrom DESeq2 estimateSizeFactors
#' @importMethodsFrom pipette droplevels2 encode export import
#' @importMethodsFrom syntactic camelCase makeNames
NULL



## S3 generics =================================================================

#' @importFrom edgeR calcNormFactors cpm scaleOffset
NULL



## Standard functions ==========================================================

#' @importFrom AcidBase dots methodFormals realpath showSlotInfo standardizeCall
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning dl h1 h2
#' h3 toInlineString txt ul
#' @importFrom AcidExperiment detectLanes importSampleData
#' makeSummarizedExperiment matchInterestingGroups
#' @importFrom AcidGenomes detectOrganism emptyRanges importTx2Gene
#' makeGRangesFromEnsembl makeGRangesFromGFF
#' @importFrom AcidPlots !! acid_coord_flip acid_geom_abline acid_geom_bar
#' acid_geom_label_repel acid_scale_y_continuous_nopad autoDiscreteColorScale
#' autoDiscreteFillScale matchLabels pretty_breaks  sym wrap_plots
#' @importFrom DESeq2 DESeq DESeqDataSet fpkm rlog
#' varianceStabilizingTransformation
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom bcbioBase getGTFFileFromYAML getMetricsFromYAML
#' getSampleDataFromYAML importDataVersions importProgramVersions projectDir
#' runDate sampleDirs
#' @importFrom edgeR DGEList
#' @importFrom ggplot2 aes expand_limits geom_point geom_smooth ggplot ggtitle
#' guides labs scale_y_continuous theme
#' @importFrom goalie areDisjointSets areIntersectingSets areSetEqual assert
#' bapply hasDimnames hasLength hasUniqueCols hasValidDimnames isADirectory
#' isAFile isAURL isAny isCharacter isDirectory isFile isFlag isGGScale isInt
#' isInRange isNonNegative isNumber isPositive isProportion isString isSubset
#' validNames validate validateClasses
#' @importFrom methods as as<- is new setAs setClass setMethod setValidity slot
#' slot<- validObject .hasSlot
#' @importFrom sessioninfo session_info
#' @importFrom tximport tximport
#' @importFrom utils capture.output packageName packageVersion
NULL
