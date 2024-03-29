#' bcbioRNASeq
#'
#' Import and analyze bcbio RNA-seq data.
#'
#' @aliases NULL
#' @keywords internal
"_PACKAGE"



## S4 classes ==================================================================

#' @importClassesFrom AcidGenomes TxToGene
#' @importClassesFrom DESeq2 DESeqDataSet DESeqTransform
#' @importClassesFrom IRanges DFrameList
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
#' plotCountsPerFeature plotDegHeatmap plotExonicMappingRate
#' plotFeatureSaturation plotGenderMarkers plotIntronicMappingRate
#' plotMappedReads plotMappingRate plotPca plotPcaCovariates plotMeanSd
#' plotPseudoVsAlignedCounts plotQc plotRrnaMappingRate plotTotalReads
#' relativeLogExpression sampleData showHeader slotAlignedCounts
#' stripTranscriptVersions tmm
#' @importFrom AcidGenomes TxToGene
#' @importFrom BiocGenerics as.data.frame counts counts<- do.call
#' estimateSizeFactors plotDispEsts updateObject width
#' @importFrom S4Vectors Rle cbind droplevels getListElement head lapply mcols
#' mcols<- metadata metadata<-
#' @importFrom SummarizedExperiment assay assay<- assayNames assayNames<- assays
#' assays<- colData colData<- rowData rowRanges
#' @importFrom methods coerce show
#' @importFrom pipette export import
NULL

#' @importMethodsFrom AcidBase showHeader
#' @importMethodsFrom AcidExperiment droplevels2 humanize interestingGroups
#' interestingGroups<- mapGenesToRownames metrics sampleData
#' stripTranscriptVersions
#' @importMethodsFrom AcidGenomes stripTranscriptVersions
#' @importMethodsFrom AcidPlots plotCorrelationHeatmap plotCountsCorrelation
#' plotCountsCorrelationHeatmap plotCountsPerBiotype plotCountsPerBroadClass
#' plotCountsPerFeature plotGenderMarkers plotPca
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
#' @importFrom AcidGenomes detectOrganism emptyRanges importTxToGene
#' makeGRangesFromEnsembl makeGRangesFromGff
#' @importFrom AcidPlots .data acid_discrete_coord_flip acid_geom_abline
#' acid_geom_bar acid_geom_label_repel acid_scale_color_discrete
#' acid_scale_fill_discrete acid_scale_y_continuous_nopad matchLabels
#' pretty_breaks wrap_plots
#' @importFrom DESeq2 DESeq DESeqDataSet fpkm rlog
#' varianceStabilizingTransformation
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom bcbioBase getGtfFileFromYaml getMetricsFromYaml
#' getSampleDataFromYaml importDataVersions importProgramVersions projectDir
#' runDate sampleDirs
#' @importFrom edgeR DGEList
#' @importFrom ggplot2 aes expand_limits geom_point geom_smooth ggplot ggtitle
#' guides labs scale_y_continuous theme
#' @importFrom goalie allAreMatchingFixed areDisjointSets areIntersectingSets
#' areSetEqual assert bapply hasDimnames hasLength hasUniqueCols
#' hasValidDimnames isADir isAFile isAUrl isAny isCharacter isDir
#' isFile isFlag isGgscale isInt isInRange isNonNegative isNumber isPositive
#' isProportion isString isSubset isTximport requireNamespaces validNames
#' validate validateClasses
#' @importFrom methods as as<- is new setAs setClass setMethod setValidity slot
#' slot<- validObject .hasSlot
#' @importFrom tximport tximport
#' @importFrom utils capture.output packageName packageVersion sessionInfo
NULL
