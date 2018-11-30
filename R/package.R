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
#' @importFrom assertive.base assert_all_are_true assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_greater_than
#'   assert_all_are_in_left_open_range assert_all_are_in_range
#'   assert_all_are_non_negative assert_all_are_positive is_positive
#' @importFrom assertive.properties assert_has_colnames assert_has_dimnames
#'   assert_has_dims assert_has_names assert_has_no_duplicates assert_has_rows
#'   assert_is_non_empty assert_is_null assert_is_of_length assert_is_scalar
#'   assert_is_vector has_dimnames has_dims has_names
#' @importFrom assertive.sets are_disjoint_sets assert_are_disjoint_sets
#'   assert_are_intersecting_sets assert_are_set_equal assert_is_subset
#'   is_subset
#' @importFrom assertive.strings assert_all_are_matching_regex
#'   assert_all_are_non_empty_character
#' @importFrom assertive.types assert_is_a_bool assert_is_a_number
#'   assert_is_a_string assert_is_all_of assert_is_an_integer assert_is_any_of
#'   assert_is_character assert_is_data.frame assert_is_factor assert_is_formula
#'   assert_is_integer assert_is_list assert_is_matrix assert_is_numeric
#'   assert_is_tbl_df is_a_string
#' @importFrom assertthat assert_that validate_that
#' @importFrom basejump Gene2Symbol Tx2Gene assertFormalGene2Symbol
#'   basejump_geom_abline basejump_geom_label basejump_geom_label_repel camel
#'   coerceS4ToList convertGenesToSymbols detectLanes detectOrganism emptyRanges
#'   formalsList import initDir interestingGroups interestingGroups<-
#'   lanePattern makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeSummarizedExperiment mapGenesToRownames markdownHeader markdownList
#'   markdownPlotlist matchArgsToDoCall matchInterestingGroups meltCounts
#'   methodFormals metrics organism plotGenesDetected plotHeatmap
#'   prepareTemplate printString realpath relevelColData relevelRowRanges
#'   removeNA sampleData sampleData<- sanitizeRowData sanitizeSampleData
#'   separator showSlotInfo snake standardizeCall stripTranscriptVersions
#'   uniteInterestingGroups validateClasses
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
#' @importFrom goalie assertAreValidNames assertHasRownames
#'   assertHasValidDimnames assertIsAlpha assertIsStringOrNULL
#'   assertIsAnImplicitInteger assertIsAnImplicitIntegerOrNULL
#'   assertIsColorScaleDiscreteOrNULL assertIsFillScaleDiscreteOrNULL
#'   assertIsHeaderLevel assertIsHexColorFunctionOrNULL assertIsImplicitInteger
#'   hasRownames hasUniqueCols
#' @importFrom grid arrow unit
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom matrixStats colMedians
#' @importFrom methods .hasSlot as as<- is new setAs show slot slot<-
#'   validObject
#' @importFrom readr read_csv read_tsv write_csv
#' @importFrom rlang !! !!! := UQ has_length sym syms
#' @importFrom scales pretty_breaks
#' @importFrom sessioninfo session_info
#' @importFrom stringr str_match str_trunc
#' @importFrom tibble as_tibble column_to_rownames remove_rownames
#'   rownames_to_column tibble
#' @importFrom tximport tximport
#' @importFrom utils capture.output globalVariables packageVersion
#' @importFrom vsn meanSdPlot
"_PACKAGE"
